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

import os
from dataclasses import dataclass

import numpy as np

COARSE_BITS = 16
COARSE_WRAP = 1 << COARSE_BITS


BLOCK_WORDS = 2_000_000      # words processed at a time (bounds temporary memory)


class WordFile:
    """The event words as a sequence of little-endian uint32, read block by
    block from disk (the whole file is never loaded).  Trailing zero padding is
    excluded.

    For a .raw file the words are the whole file.  For a streamed .obf the
    words are the event-stream stack inside the container (byte_offset,
    n_words; see obf.find_event_stream)."""

    def __init__(self, path: str, block: int | None = None,
                 byte_offset: int = 0, n_words: int | None = None):
        self.path = path
        self.block = block or BLOCK_WORDS
        self.byte_offset = byte_offset
        n_total = os.path.getsize(path) // 4 if n_words is None else n_words
        if n_total == 0:
            raise ValueError(f"{path}: file is empty")
        last = 0
        for end in range(n_total, 0, -self.block):
            s0 = max(0, end - self.block)
            nz = np.flatnonzero(self.read(s0, end))
            if nz.size:
                last = s0 + int(nz[-1]) + 1
                break
        if last == 0:
            raise ValueError(f"{path}: file contains only zeros")
        self.n_words = last
        self.n_padding = n_total - last
        self.blocks = [(s0, min(last, s0 + self.block)) for s0 in range(0, last, self.block)]

    def read(self, s0: int, s1: int) -> np.ndarray:
        return np.fromfile(self.path, dtype="<u4", count=s1 - s0,
                           offset=self.byte_offset + 4 * s0)


@dataclass
class RawEvents:
    """Decoded content of one .raw file. All coarse times are unwrapped (int64)."""
    n_words: int
    n_padding_words: int
    # photons (all elements except the sync element)
    ph_element: np.ndarray   # uint8
    ph_code: np.ndarray      # uint16, raw 10-bit TDC code
    ph_coarse: np.ndarray    # int64
    # sync records (None when decoded with keep_sync=False; use `sync` instead)
    sy_code: np.ndarray | None
    sy_coarse: np.ndarray | None
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
    n_sync: int
    sync: "SyncReader"       # block-wise access to the sync records

    @property
    def n_photons(self):
        return int(self.ph_element.size)


def is_obf(path: str) -> bool:
    with open(path, "rb") as f:
        return f.read(10) == b"OMAS_BF\n\xff\xff"


def open_words(path: str, block: int | None = None) -> WordFile:
    """WordFile for a .raw file or for the event-stream stack of a streamed .obf
    (recognised by its file signature, not by the extension)."""
    if is_obf(path):
        from .obf import find_event_stream
        es = find_event_stream(path)
        return WordFile(path, block, byte_offset=es.byte_offset, n_words=es.n_words)
    return WordFile(path, block)


def read_words(path: str) -> tuple[np.ndarray, int]:
    """Read the whole file as little-endian uint32 and strip trailing zero padding."""
    wf = open_words(path)
    return wf.read(0, wf.n_words), wf.n_padding


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


class SyncReader:
    """Block-wise access to the sync records of a .raw file, re-read from disk
    on demand.  The sync channel is ~3/4 of all words; keeping it out of memory
    is the main memory saving.  `carry[i]` is the unwrap state before word
    block i, so every block can be decoded on its own with identical results.
    """

    def __init__(self, wf: WordFile, carry: list, counts: np.ndarray, sync_element: int,
                 code_hist: np.ndarray):
        self.wf = wf
        self.carry = carry
        self.cum = np.concatenate([[0], np.cumsum(counts)]).astype(np.int64)
        self.n = int(self.cum[-1])
        self.sync_element = sync_element
        self.code_hist = code_hist          # histogram of the 10-bit TDC codes (for calibration)
        self._cache = {}

    def _block(self, i):
        if i not in self._cache:
            s0, s1 = self.wf.blocks[i]
            w = self.wf.read(s0, s1)
            coarse, _ = unwrap_coarse(w, prev=self.carry[i])
            m = ((w >> 31) == 1) & (((w >> 26) & 0x1F) == self.sync_element)
            if len(self._cache) >= 2:
                self._cache.pop(next(iter(self._cache)))
            self._cache[i] = (coarse[m], ((w[m] >> 16) & 0x3FF).astype(np.uint16))
        return self._cache[i]

    def read(self, a: int, b: int):
        """(coarse int64, code uint16) of sync records a..b-1."""
        if not 0 <= a <= b <= self.n:
            raise IndexError(f"sync records {a}..{b} outside 0..{self.n}")
        i0 = int(np.searchsorted(self.cum, a, side="right")) - 1
        parts_c, parts_k = [], []
        i = i0
        while i < len(self.wf.blocks) and self.cum[i] < b:
            c, k = self._block(i)
            lo = max(a, self.cum[i]) - self.cum[i]
            hi = min(b, self.cum[i + 1]) - self.cum[i]
            parts_c.append(c[lo:hi])
            parts_k.append(k[lo:hi])
            i += 1
        if not parts_c:
            return np.empty(0, np.int64), np.empty(0, np.uint16)
        if len(parts_c) == 1:
            return parts_c[0], parts_k[0]
        return np.concatenate(parts_c), np.concatenate(parts_k)

    def iter_blocks(self):
        """Yield (coarse, code) per word block, in stream order (not cached)."""
        for i in range(len(self.wf.blocks)):
            if self.cum[i + 1] > self.cum[i]:
                if i in self._cache:
                    yield self._cache[i]
                    continue
                s0, s1 = self.wf.blocks[i]
                w = self.wf.read(s0, s1)
                coarse, _ = unwrap_coarse(w, prev=self.carry[i])
                m = ((w >> 31) == 1) & (((w >> 26) & 0x1F) == self.sync_element)
                yield coarse[m], ((w[m] >> 16) & 0x3FF).astype(np.uint16)

    def clear_cache(self):
        self._cache.clear()


def decode(path: str, sync_element: int = 23, keep_sync: bool = True,
           block: int | None = None) -> RawEvents:
    """Decode a .raw file (or the event stream of a streamed .obf), reading it
    block by block (BLOCK_WORDS at a time).

    keep_sync=False does not keep the sync records in memory (sy_coarse and
    sy_code are None); they stay available block-wise through `ev.sync`.
    """
    wf = open_words(path, block)
    blocks = wf.blocks

    # pass 1: unwrap state per block, record counts, sync code histogram
    n_ph = n_mk = n_gen = 0
    sy_counts = np.zeros(len(blocks), np.int64)
    sy_hist = np.zeros(1024, np.int64)
    carry, prev, back = [], None, 0
    for i, (s0, s1) in enumerate(blocks):
        w = wf.read(s0, s1)
        carry.append(prev)
        coarse, b = unwrap_coarse(w, prev=prev)
        back = max(back, b)
        prev = (int(w[-1] & 0xFFFF), int(coarse[-1]))
        del coarse
        _, is_ph, is_sync, is_mk, is_gen = _classify(w, sync_element)
        n_ph += int(is_ph.sum())
        sy_counts[i] = int(is_sync.sum())
        sy_hist += np.bincount((w[is_sync] >> 16) & 0x3FF, minlength=1024)
        n_mk += int(is_mk.sum())
        n_gen += int(is_gen.sum())
    n_sy = int(sy_counts.sum())
    sync = SyncReader(wf, carry, sy_counts, sync_element, sy_hist)

    ph_element = np.empty(n_ph, np.uint8)
    ph_code = np.empty(n_ph, np.uint16)
    ph_coarse = np.empty(n_ph, np.int64)
    sy_code = np.empty(n_sy, np.uint16) if keep_sync else None
    sy_coarse = np.empty(n_sy, np.int64) if keep_sync else None
    mk_w = np.empty(n_mk, np.uint32)
    mk_coarse = np.empty(n_mk, np.int64)
    mk_word = np.empty(n_mk, np.int64)
    gen_w = np.empty(n_gen, np.uint32)
    gen_coarse = np.empty(n_gen, np.int64)

    # pass 2: fill
    i_ph = i_sy = i_mk = i_gen = 0
    for i, (s0, s1) in enumerate(blocks):
        w = wf.read(s0, s1)
        coarse, _ = unwrap_coarse(w, prev=carry[i])
        element, is_ph, is_sync, is_mk, is_gen = _classify(w, sync_element)
        k = int(is_ph.sum())
        ph_element[i_ph:i_ph + k] = element[is_ph]
        ph_code[i_ph:i_ph + k] = (w[is_ph] >> 16) & 0x3FF
        ph_coarse[i_ph:i_ph + k] = coarse[is_ph]
        i_ph += k
        if keep_sync:
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
        del w, coarse, element, is_ph, is_sync, is_mk, is_gen

    return RawEvents(
        n_words=int(wf.n_words), n_padding_words=wf.n_padding,
        ph_element=ph_element, ph_code=ph_code, ph_coarse=ph_coarse,
        sy_code=sy_code, sy_coarse=sy_coarse,
        mk_bits=((mk_w >> 26) & 0x7).astype(np.uint8),
        mk_payload=((mk_w >> 16) & 0x3FF).astype(np.uint16), mk_coarse=mk_coarse, mk_word=mk_word,
        gen_type=((gen_w >> 26) & 0x3F).astype(np.uint8),
        gen_payload=((gen_w >> 16) & 0x3FF).astype(np.uint16), gen_coarse=gen_coarse,
        max_backstep=back, sync_element=sync_element, n_sync=n_sy, sync=sync,
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
        "sync_records": int(ev.n_sync),
        "markers_by_bits(frame,line,pixel)": {f"{int(a):03b}": int(b) for a, b in zip(mb, mc)},
        "general_messages_by_type": {int(a): int(b) for a, b in zip(gt, gc)},
        "max_coarse_backstep_ticks": ev.max_backstep,
    }
