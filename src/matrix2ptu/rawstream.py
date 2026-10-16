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


BLOCK_WORDS = 8_000_000      # words processed at a time (bounds temporary memory)


@dataclass
class RawEvents:
    """Decoded content of one .raw file. All coarse times are unwrapped (int64)."""
    n_words: int
    n_padding_words: int
    # photons (all elements except the sync element)
    ph_element: np.ndarray   # uint8
    ph_code: np.ndarray      # uint16, raw 10-bit TDC code
    ph_coarse: np.ndarray    # int64
    # sync records
    sy_code: np.ndarray
    sy_coarse: np.ndarray
    # scan markers
    mk_bits: np.ndarray      # uint8: bit2=frame, bit1=line, bit0=pixel
    mk_payload: np.ndarray
    mk_coarse: np.ndarray
    mk_word: np.ndarray      # int64, word index in file (for traceability)
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
    # search for the last non-zero word from the end, in blocks (no full-size index array)
    last = 0
    for end in range(w.size, 0, -BLOCK_WORDS):
        nz = np.flatnonzero(w[max(0, end - BLOCK_WORDS):end])
        if nz.size:
            last = max(0, end - BLOCK_WORDS) + int(nz[-1]) + 1
            break
    if last == 0:
        raise ValueError(f"{path}: file contains only zeros")
    return w[:last], int(w.size - last)


def unwrap_coarse(words: np.ndarray, max_back: int = 1024,
                  prev: tuple[int, int] | None = None) -> tuple[np.ndarray, int]:
    """Unwrap the 16-bit coarse counter of all words into an int64 time axis.

    The difference of consecutive 16-bit values, taken modulo 2**16 and read
    as a signed int16, is the true step as long as no gap between consecutive
    words exceeds 32767 ticks (~0.33 ms).  The dense sync channel (a record
    every ~78 ticks) guarantees this; build_pulse_train() additionally checks
    the sync spacing.

    For block-wise processing, `prev` = (16-bit value, unwrapped value) of the
    word preceding this block; the result is identical to unwrapping the whole
    stream at once.
    """
    c16 = (words & 0xFFFF).astype(np.uint16)
    if prev is None:
        d = np.diff(c16).view(np.int16)          # modular difference, signed
        base = int(c16[0])
    else:
        d = np.diff(c16, prepend=np.uint16(prev[0])).view(np.int16)
        base = prev[1]
    back = int(-d.min()) if d.size and d.min() < 0 else 0
    if back > max_back:
        raise ValueError(
            f"coarse counter steps back by {back} ticks between consecutive words; "
            "the stream does not look time-ordered - cannot unwrap safely")
    out = np.empty(words.size, dtype=np.int64)
    if prev is None:
        out[0] = base
        np.cumsum(d, dtype=np.int64, out=out[1:])
        out[1:] += base
    else:
        np.cumsum(d, dtype=np.int64, out=out)
        out += base
    return out, back


def _classify(w: np.ndarray, sync_element: int):
    hi = (w >> 26).astype(np.uint8)          # top 6 bits = message type
    is_ph = hi >= 32
    element = hi & 0x1F
    is_sync = is_ph & (element == sync_element)
    is_ph &= ~is_sync
    is_mk = (hi >> 3) == 1                   # 0b001xxx
    is_gen = hi < 8                          # 0b000xxx
    return element, is_ph, is_sync, is_mk, is_gen


def decode(path: str, sync_element: int = 23) -> RawEvents:
    """Decode a .raw file.  Works block-wise (BLOCK_WORDS at a time) so that
    temporary arrays stay small; the output arrays are allocated once."""
    words, npad = read_words(path)
    n = words.size
    blocks = [(s0, min(n, s0 + BLOCK_WORDS)) for s0 in range(0, n, BLOCK_WORDS)]

    # pass 1: count record types to allocate the outputs exactly
    n_ph = n_sy = n_mk = n_gen = 0
    for s0, s1 in blocks:
        _, is_ph, is_sync, is_mk, is_gen = _classify(words[s0:s1], sync_element)
        n_ph += int(is_ph.sum())
        n_sy += int(is_sync.sum())
        n_mk += int(is_mk.sum())
        n_gen += int(is_gen.sum())
    ph_element = np.empty(n_ph, np.uint8)
    ph_code = np.empty(n_ph, np.uint16)
    ph_coarse = np.empty(n_ph, np.int64)
    sy_code = np.empty(n_sy, np.uint16)
    sy_coarse = np.empty(n_sy, np.int64)
    mk_w = np.empty(n_mk, np.uint32)
    mk_coarse = np.empty(n_mk, np.int64)
    mk_word = np.empty(n_mk, np.int64)
    gen_w = np.empty(n_gen, np.uint32)
    gen_coarse = np.empty(n_gen, np.int64)

    # pass 2: unwrap the coarse counter (carried across blocks) and fill
    i_ph = i_sy = i_mk = i_gen = 0
    back = 0
    prev = None
    for s0, s1 in blocks:
        w = words[s0:s1]
        coarse, b = unwrap_coarse(w, prev=prev)
        back = max(back, b)
        prev = (int(w[-1] & 0xFFFF), int(coarse[-1]))
        element, is_ph, is_sync, is_mk, is_gen = _classify(w, sync_element)
        k = int(is_ph.sum())
        ph_element[i_ph:i_ph + k] = element[is_ph]
        ph_code[i_ph:i_ph + k] = (w[is_ph] >> 16) & 0x3FF
        ph_coarse[i_ph:i_ph + k] = coarse[is_ph]
        i_ph += k
        k = int(is_sync.sum())
        sy_code[i_sy:i_sy + k] = (w[is_sync] >> 16) & 0x3FF
        sy_coarse[i_sy:i_sy + k] = coarse[is_sync]
        i_sy += k
        k = int(is_mk.sum())
        mk_w[i_mk:i_mk + k] = w[is_mk]
        mk_coarse[i_mk:i_mk + k] = coarse[is_mk]
        mk_word[i_mk:i_mk + k] = np.flatnonzero(is_mk) + s0
        i_mk += k
        k = int(is_gen.sum())
        gen_w[i_gen:i_gen + k] = w[is_gen]
        gen_coarse[i_gen:i_gen + k] = coarse[is_gen]
        i_gen += k
        del coarse, element, is_ph, is_sync, is_mk, is_gen

    return RawEvents(
        n_words=int(n), n_padding_words=npad,
        ph_element=ph_element, ph_code=ph_code, ph_coarse=ph_coarse,
        sy_code=sy_code, sy_coarse=sy_coarse,
        mk_bits=((mk_w >> 26) & 0x7).astype(np.uint8),
        mk_payload=((mk_w >> 16) & 0x3FF).astype(np.uint16), mk_coarse=mk_coarse, mk_word=mk_word,
        gen_type=((gen_w >> 26) & 0x3F).astype(np.uint8),
        gen_payload=((gen_w >> 16) & 0x3FF).astype(np.uint16), gen_coarse=gen_coarse,
        max_backstep=back, sync_element=sync_element,
    )


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
