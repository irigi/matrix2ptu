"""Decoding of the Abberior MATRIX ".raw" stream (32-bit messages, "V5").

Word layout (from Abberior's description, verified against the data where
noted in docs/REPORT.md):

    bit 31 = 1      photon / sync     [1:1][5:element][10:tdc][16:coarse]
    bits 31..29=001 scan marker       [3:001][1:frame][1:line][1:pixel][10:payload][16:coarse]
    bits 31..29=000 "general" message [6:type][10:payload][16:coarse]

Element numbers 0..22 are the 23 SPAD pixels; element 23 is the laser-sync
channel (one record per 32 laser pulses).  These assignments were established
from the data (counts/periodicity), not from a written specification, and can be
overridden via the configuration.

The 16-bit coarse counter wraps without explicit overflow records.  Because the
sync channel produces a record every ~78 coarse ticks, consecutive words are
never more than a few ticks apart, so the counter can be unwrapped by taking the
signed modular difference between consecutive words.  The stream is almost
time-ordered: steps of -1/-2 ticks occur (records from different TDCs are
merged with small reordering).  The unwrapping tolerates this.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

COARSE_BITS = 16
COARSE_WRAP = 1 << COARSE_BITS


@dataclass
class RawEvents:
    """Decoded content of one .raw file. All coarse times are unwrapped (int64)."""
    n_words: int
    n_padding_words: int
    # photons (all elements except the sync element)
    ph_element: np.ndarray   # uint8
    ph_code: np.ndarray      # uint16, raw 10-bit TDC code
    ph_coarse: np.ndarray    # int64
    ph_word: np.ndarray      # int64, word index in file (for traceability)
    # sync records
    sy_code: np.ndarray
    sy_coarse: np.ndarray
    # scan markers
    mk_bits: np.ndarray      # uint8: bit2=frame, bit1=line, bit0=pixel
    mk_payload: np.ndarray
    mk_coarse: np.ndarray
    mk_word: np.ndarray
    # general messages (type 0..7), kept for reporting
    gen_type: np.ndarray
    gen_payload: np.ndarray
    gen_coarse: np.ndarray
    max_backstep: int        # largest backward coarse step seen (ticks)
    sync_element: int

    @property
    def n_photons(self):
        return int(self.ph_element.size)


def read_words(path: str) -> tuple[np.ndarray, int]:
    """Read the file as little-endian uint32 and strip trailing zero padding."""
    w = np.fromfile(path, dtype="<u4")
    if w.size == 0:
        raise ValueError(f"{path}: file is empty")
    nz = np.flatnonzero(w)
    if nz.size == 0:
        raise ValueError(f"{path}: file contains only zeros")
    last = nz[-1] + 1
    return w[:last], int(w.size - last)


def unwrap_coarse(words: np.ndarray, max_back: int = 1024) -> tuple[np.ndarray, int]:
    """Unwrap the 16-bit coarse counter of all words into an int64 time axis.

    The difference of consecutive 16-bit values, taken modulo 2**16 and read
    as a signed int16, is the true step as long as no gap between consecutive
    words exceeds 32767 ticks (~0.33 ms).  The dense sync channel (a record
    every ~78 ticks) guarantees this; build_pulse_train() additionally checks
    the sync spacing.
    """
    c16 = (words & 0xFFFF).astype(np.uint16)
    d = np.diff(c16).view(np.int16)          # modular difference, signed
    back = int(-d.min()) if d.size and d.min() < 0 else 0
    if back > max_back:
        raise ValueError(
            f"coarse counter steps back by {back} ticks between consecutive words; "
            "the stream does not look time-ordered - cannot unwrap safely")
    out = np.empty(words.size, dtype=np.int64)
    out[0] = c16[0]
    np.cumsum(d, dtype=np.int64, out=out[1:])
    out[1:] += int(c16[0])
    return out, back


def decode(path: str, sync_element: int = 23) -> RawEvents:
    words, npad = read_words(path)
    coarse, back = unwrap_coarse(words)

    hi = (words >> 26).astype(np.uint8)        # top 6 bits = message type
    is_ph = hi >= 32
    element = hi & 0x1F
    is_sync = is_ph & (element == sync_element)
    idx_ph = np.flatnonzero(is_ph & ~is_sync)
    idx_sy = np.flatnonzero(is_sync)
    del is_ph, is_sync
    idx_mk = np.flatnonzero((hi >> 3) == 1)    # 0b001xxx
    idx_gen = np.flatnonzero(hi < 8)           # 0b000xxx
    del hi
    element = element[idx_ph]

    def code(idx):
        return ((words[idx] >> 16) & 0x3FF).astype(np.uint16)
    ev = RawEvents(
        n_words=int(words.size), n_padding_words=npad,
        ph_element=element, ph_code=code(idx_ph), ph_coarse=coarse[idx_ph],
        ph_word=idx_ph.astype(np.int64),
        sy_code=code(idx_sy), sy_coarse=coarse[idx_sy],
        mk_bits=((words[idx_mk] >> 26) & 0x7).astype(np.uint8),
        mk_payload=code(idx_mk), mk_coarse=coarse[idx_mk], mk_word=idx_mk.astype(np.int64),
        gen_type=((words[idx_gen] >> 26) & 0x3F).astype(np.uint8),
        gen_payload=code(idx_gen), gen_coarse=coarse[idx_gen],
        max_backstep=back, sync_element=sync_element,
    )
    return ev


def summarize(ev: RawEvents) -> dict:
    el, cnt = np.unique(ev.ph_element, return_counts=True)
    mb, mc = np.unique(ev.mk_bits, return_counts=True)
    gt, gc = np.unique(ev.gen_type, return_counts=True)
    return {
        "words": ev.n_words,
        "padding_words": ev.n_padding_words,
        "photons": ev.n_photons,
        "photons_per_element": {int(a): int(b) for a, b in zip(el, cnt)},
        "sync_records": int(ev.sy_coarse.size),
        "markers_by_bits(frame,line,pixel)": {f"{int(a):03b}": int(b) for a, b in zip(mb, mc)},
        "general_messages_by_type": {int(a): int(b) for a, b in zip(gt, gc)},
        "max_coarse_backstep_ticks": ev.max_backstep,
    }
