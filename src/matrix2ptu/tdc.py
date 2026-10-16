"""TDC fine-time calibration by the statistical code-density method.

What the data show (docs/REPORT.md, "TDC code structure"):
  * The 10-bit TDC code is [2-bit quarter][8-bit delay-line tap].
  * The quarter q = code >> 8 adds q/4 of a coarse tick.
  * Inside a quarter the tap number *decreases* as the event time increases
    (the delay line measures time to the next clock phase).
  * Only ~145 of the 256 taps are used and their widths differ strongly
    (large differential non-linearity), so a linear code->time mapping is wrong.

Code-density calibration: if events are uniformly distributed in phase relative
to the coarse clock (true here: laser and coarse clock are not frequency-locked,
ratio 78.334/32 ticks per pulse), the number of events in a code is proportional
to that code's width.  The cumulative histogram therefore maps each code to its
position inside the quarter.  Each element has its own TDC, so each element
(and the sync channel) is calibrated separately.

This replaces Abberior's internal (unpublished) per-TDC calibration.  Its
accuracy is checked in validation via sync periodicity and by comparing decay
histograms with those Imspector stored in the .msr file.
"""
from __future__ import annotations

import numpy as np

N_CODES = 1024
TAPS = 256
QUARTERS = 4


def code_density_table(codes: np.ndarray, reverse_taps: bool = True,
                       min_counts: int = 10000) -> np.ndarray:
    """Return a 1024-entry table: code -> fractional coarse tick in [0, 1).

    The value is the centre of the code's bin.  Unused codes get NaN.
    """
    h = np.bincount(codes.astype(np.int64), minlength=N_CODES).astype(float)
    if h.sum() < min_counts:
        raise ValueError(
            f"only {int(h.sum())} events available for TDC calibration (need >= {min_counts}); "
            "supply a calibration table from a larger dataset")
    table = np.full(N_CODES, np.nan)
    for q in range(QUARTERS):
        hq = h[q * TAPS:(q + 1) * TAPS]
        tot = hq.sum()
        if tot == 0:
            continue
        if reverse_taps:
            r = hq[::-1]
            pos = ((np.cumsum(r) - r / 2) / tot)[::-1]
        else:
            pos = (np.cumsum(hq) - hq / 2) / tot
        pos[hq == 0] = np.nan
        table[q * TAPS:(q + 1) * TAPS] = (q + pos) / QUARTERS
    return table


class TdcCalibration:
    """Per-channel code -> fraction-of-coarse-tick tables."""

    def __init__(self, tables: dict[int, np.ndarray], method: str):
        self.tables = tables
        self.method = method

    @classmethod
    def from_events(cls, ph_element, ph_code, sy_code, sync_element=23,
                    reverse_taps=True):
        tables = {}
        for e in np.unique(ph_element):
            tables[int(e)] = code_density_table(ph_code[ph_element == e], reverse_taps,
                                                min_counts=1000)
        tables[sync_element] = code_density_table(sy_code, reverse_taps)
        return cls(tables, "code-density (self-calibrated from this file)")

    def fine_one(self, element: int, code: np.ndarray) -> np.ndarray:
        """Fine time for codes that all come from one element (no masks)."""
        out = self.tables[int(element)][code]
        if np.isnan(out).any():
            raise ValueError(f"element {element}: TDC codes without calibration value")
        return out

    def fine(self, element: np.ndarray, code: np.ndarray) -> np.ndarray:
        out = np.empty(code.size)
        for e, tab in self.tables.items():
            m = element == e
            if m.any():
                out[m] = tab[code[m]]
        if np.isnan(out).any():
            raise ValueError("events with TDC codes that have no calibration value")
        return out

    def save(self, path):
        np.savez(path, method=self.method,
                 **{f"el{e:02d}": t for e, t in self.tables.items()})

    @classmethod
    def load(cls, path):
        z = np.load(path)
        tables = {int(k[2:]): z[k] for k in z.files if k.startswith("el")}
        return cls(tables, str(z["method"]) + f" (loaded from {path})")
