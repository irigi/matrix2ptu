"""Pixel reassignment of MATRIX detector-resolved PTU data (event level).

Works on the per-sequence detector-resolved PTU files written by the converter
(``split_sequences: true``).  Each photon record keeps its detector element
(PTU channel), its micro time (dtime) and its macro time (nsync); the image
position of a photon follows from the PTU line markers:

    line pass k = last LineStart before the photon
    x  = (t - LineStart_k) / (LineStop_k - LineStart_k) * pixels_x   (continuous)
    y  = k // line_accumulation,   pass within the line = k % line_accumulation

Coordinates used throughout: x = pixel column (scan direction, right),
y = line index (down), both in pixels of the acquisition grid; origin at the
top-left corner of the first pixel.  A *shift* s_e = (sx, sy) of element e
means: the element image is the reference image displaced by +s_e,
I_e(r) ~ I_ref(r - s_e).  Reassignment moves every photon of element e by
-s_e.

Reassignment inside a PTU is done by rewriting the macro time only:
    x' = x - sx           ->  new time inside the target line pass (exact to
                              1 laser period = 1/4000 px at 100 us dwell)
    y' = round(y - sy)    ->  same pass index in the target line
The micro time and the channel of every photon are copied unchanged.  The
macro time of a reassigned photon is therefore a *position code*, not its
detection time; the line/frame markers keep their original times.
"""
from __future__ import annotations

import json
import os
import struct
from dataclasses import dataclass, field

import numpy as np

from .validate import iter_ptu_records

N_EL = 23


# --------------------------------------------------------------------------- photons
@dataclass
class PhotonData:
    path: str
    pixels_x: int
    pixels_y: int
    accumulation: int
    pixel_size_nm: float
    t: np.ndarray            # nsync (int64), all photons, file order
    ch: np.ndarray           # element (uint8)
    dtime: np.ndarray        # micro time bin (uint16)
    ls: np.ndarray           # LineStart nsync per line pass
    le: np.ndarray           # LineStop nsync per line pass
    frame: np.ndarray        # Frame marker nsync
    k: np.ndarray = field(init=False)       # line pass index (-1 before the first line)
    inside: np.ndarray = field(init=False)  # photon inside a line window
    x: np.ndarray = field(init=False)       # continuous x (pixels), valid where inside
    y: np.ndarray = field(init=False)       # line index, valid where inside

    def __post_init__(self):
        k = np.searchsorted(self.ls, self.t, side="right") - 1
        kk = np.clip(k, 0, None)
        inside = (k >= 0) & (self.t < self.le[kk])
        self.k = k
        self.x = (self.t - self.ls[kk]) / (self.le[kk] - self.ls[kk]) * self.pixels_x
        self.y = kk // self.accumulation
        inside &= self.y < self.pixels_y
        self.inside = inside

    def images(self, x=None, y=None, sel=None, dtype=np.float64):
        """Per-element images (23, Y, X) from (optionally modified) positions."""
        x = self.x if x is None else x
        y = self.y if y is None else y
        m = self.inside if sel is None else (self.inside & sel)
        xi = np.floor(x[m]).astype(np.int64)
        yi = np.asarray(y[m]).astype(np.int64)
        ok = (xi >= 0) & (xi < self.pixels_x) & (yi >= 0) & (yi < self.pixels_y)
        idx = (self.ch[m][ok].astype(np.int64) * self.pixels_y + yi[ok]) * self.pixels_x + xi[ok]
        n = N_EL * self.pixels_y * self.pixels_x
        return np.bincount(idx, minlength=n).reshape(N_EL, self.pixels_y, self.pixels_x).astype(dtype)


def load_photons(ptu_path: str, accumulation: int, pixels_x: int, pixels_y: int,
                 pixel_size_nm: float) -> PhotonData:
    ts, cs, ds, ls, le, fr = [], [], [], [], [], []
    for _, d, t in iter_ptu_records(ptu_path):
        ph = (d["channel"] >= 0) & (d["marker"] == 0)
        ts.append(t[ph]); cs.append(d["channel"][ph].astype(np.uint8))
        ds.append(d["dtime"][ph].astype(np.uint16))
        ls.append(t[(d["marker"] & 1) > 0]); le.append(t[(d["marker"] & 2) > 0])
        fr.append(t[(d["marker"] & 4) > 0])
    c = np.concatenate
    return PhotonData(ptu_path, pixels_x, pixels_y, accumulation, pixel_size_nm,
                      c(ts), c(cs), c(ds), c(ls), c(le), c(fr))


# --------------------------------------------------------------------------- reassignment
def reassign_positions(pd: PhotonData, shifts_px: np.ndarray):
    """New continuous positions x' = x - sx[e], y' = y - sy[e] (no rounding)."""
    e = pd.ch.astype(np.int64)
    return pd.x - shifts_px[e, 0], pd.y - shifts_px[e, 1]


def reassign_times(pd: PhotonData, shifts_px: np.ndarray, y_rounding: str = "nearest", seed: int = 0):
    """Event-level reassignment.  Returns (new nsync, keep mask, stats, applied y shift per photon).

    Photons inside a line window are moved to x' = x - sx (exact up to one laser
    period) and to a whole line at the same pass index:
      y_rounding='nearest':    line y - round(sy)  (deterministic; residual
                               |sy - round(sy)| <= 0.5 px per element)
      y_rounding='stochastic': line y - floor(sy) - B, B ~ Bernoulli(frac(sy))
                               (unbiased mean shift; adds a 2-line spread of
                               variance frac(1-frac) px^2; fixed seed)
    Photons that would leave the field are dropped (keep = False).  Photons
    outside all line windows (flyback) keep their time: they have no image
    position before or after."""
    e = pd.ch.astype(np.int64)
    sy = shifts_px[:, 1]
    if y_rounding == "nearest":
        sy_int = np.rint(sy).astype(np.int64)
        dy = sy_int[e]
    elif y_rounding == "stochastic":
        fl = np.floor(sy).astype(np.int64)
        frac = sy - fl
        rng = np.random.default_rng(seed)
        dy = fl[e] + (rng.random(e.size) < frac[e])
        sy_int = np.rint(sy).astype(np.int64)
    else:
        raise ValueError(f"y_rounding must be 'nearest' or 'stochastic', not {y_rounding!r}")
    xn = pd.x - shifts_px[e, 0]
    yn = pd.y - dy
    acc = pd.accumulation
    kn = yn * acc + np.clip(pd.k, 0, None) % acc
    out = (xn < 0) | (xn >= pd.pixels_x) | (yn < 0) | (yn >= pd.pixels_y)
    move = pd.inside & ~out
    kn_c = np.clip(kn, 0, pd.ls.size - 1)
    span = (pd.le[kn_c] - pd.ls[kn_c]).astype(np.float64)
    tn = pd.t.copy()
    tn[move] = pd.ls[kn_c[move]] + np.floor(xn[move] / pd.pixels_x * span[move]).astype(np.int64)
    keep = ~(pd.inside & out)
    stats = {
        "photons_total": int(pd.t.size),
        "photons_in_image": int(pd.inside.sum()),
        "photons_moved": int(move.sum()),
        "photons_dropped_out_of_field": int((pd.inside & out).sum()),
        "photons_outside_line_windows_unchanged": int((~pd.inside).sum()),
        "dropped_per_element": np.bincount(e[pd.inside & out], minlength=N_EL).tolist(),
        "y_rounding": y_rounding,
        "y_shift_applied_lines": sy_int.tolist() if y_rounding == "nearest" else
        (np.bincount(e[move], weights=dy[move], minlength=N_EL) /
         np.maximum(np.bincount(e[move], minlength=N_EL), 1)).round(4).tolist(),
        "y_rounding_residual_px": (sy - sy_int).round(4).tolist() if y_rounding == "nearest" else "0 on average",
        "x_quantum_px": float(pd.pixels_x / np.median(pd.le - pd.ls)),
    }
    return tn, keep, stats, dy


# --------------------------------------------------------------------------- PTU header splice
_TY_EMPTY, _TY_STRINGS = 0xFFFF0008, (0x4001FFFF, 0x4002FFFF, 0xFFFFFFFF)


def _walk_tags(hdr: bytes):
    """Yield (ident, idx, typ, start, end, value_offset) over a PTU header."""
    pos = 16
    while True:
        ident = hdr[pos:pos + 32].split(b"\0")[0].decode()
        idx, typ = struct.unpack_from("<iI", hdr, pos + 32)
        val_off = pos + 40
        end = pos + 48
        if typ in _TY_STRINGS:
            end += struct.unpack_from("<q", hdr, val_off)[0]
        yield ident, idx, typ, pos, end, val_off
        pos = end
        if ident == "Header_End":
            return


def read_header_bytes(path: str) -> bytes:
    with open(path, "rb") as f:
        blob = f.read(1 << 20)
    for ident, _, _, _, end, _ in _walk_tags(blob):
        if ident == "Header_End":
            return blob[:end]
    raise ValueError(f"{path}: Header_End not found")


def write_derived_ptu(src_ptu: str, dst_ptu: str, extra_tags: list[tuple],
                      nsync, channel, dtime, special, replace_comment: str | None = None) -> int:
    """Write a PTU whose header is the source header verbatim plus extra tags
    (inserted before Header_End), with the record count patched and the
    File_GUID renewed.  Events must be time-sorted."""
    from .ptu_writer import _tag, encode_t3_fast, new_guid

    hdr = read_header_bytes(src_ptu)
    parts = [hdr[:16]]
    replaced = {(t[0], t[2] if len(t) > 2 else -1) for t in extra_tags}
    for ident, idx, typ, start, end, _ in _walk_tags(hdr):
        if (ident, idx) in replaced:
            continue                      # superseded by an extra tag
        if ident == "Header_End":
            parts.extend(_tag(*t) if len(t) > 2 else _tag(t[0], t[1]) for t in extra_tags)
            parts.append(hdr[start:end])
        elif ident == "File_GUID":
            parts.append(_tag("File_GUID", new_guid()))
        elif ident == "File_Comment" and replace_comment is not None:
            parts.append(_tag("File_Comment", replace_comment))
        else:
            parts.append(hdr[start:end])
    new_hdr = b"".join(parts)
    nrec_off = next(v for i, _, _, _, _, v in _walk_tags(new_hdr) if i == "TTResult_NumberOfRecords")
    nrec, hi_prev = 0, 0
    with open(dst_ptu, "wb") as f:
        f.write(new_hdr)
        for s0 in range(0, len(nsync), 4_000_000):
            sl = slice(s0, s0 + 4_000_000)
            rec = encode_t3_fast(nsync[sl], channel[sl], dtime[sl], special[sl], hi_prev)
            rec.astype("<u4", copy=False).tofile(f)
            nrec += rec.size
            hi_prev = int(nsync[sl][-1]) // 1024
        f.seek(nrec_off)
        f.write(struct.pack("<q", nrec))
    return nrec


def merge_events(pd: PhotonData, t_ph, keep, pooled: bool):
    """Time-sort photons (new times) with the original markers; markers first at
    equal time (as the converter does)."""
    from .converter import MARK_FRAME, MARK_LINE_START, MARK_LINE_STOP
    tp = t_ph[keep]
    cp = np.zeros(tp.size, np.uint8) if pooled else pd.ch[keep]
    dp = pd.dtime[keep]
    mk_t = np.concatenate([pd.frame, pd.ls, pd.le]).astype(np.int64)
    mk_c = np.concatenate([np.full(pd.frame.size, MARK_FRAME), np.full(pd.ls.size, MARK_LINE_START),
                           np.full(pd.le.size, MARK_LINE_STOP)]).astype(np.uint8)
    n = np.concatenate([mk_t, tp])
    c = np.concatenate([mk_c, cp])
    d = np.concatenate([np.zeros(mk_t.size, np.uint16), dp])
    s = np.concatenate([np.ones(mk_t.size, np.uint8), np.zeros(tp.size, np.uint8)])
    o = np.lexsort((1 - s, n))          # by time, markers before photons
    return n[o], c[o], d[o], s[o]


def save_json(obj, path):
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w") as f:
        json.dump(obj, f, indent=2, default=lambda v: v.tolist() if hasattr(v, "tolist") else str(v))
