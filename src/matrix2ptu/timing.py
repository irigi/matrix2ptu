"""Reconstruction of the laser pulse train and photon micro/macro times.

The MATRIX does not deliver PicoQuant-style (nsync, dtime) pairs.  It delivers
time stamps on its own coarse clock (~97.92 MHz, measured) plus a TDC fine
time, and a "sync" record every 32nd laser pulse (Abberior: "the sync channel
is just divided by 32 and tagged like a photon channel").

Procedure (Abberior describes the same principle: "photon times relative to
the sync pulse are calculated by tracking phase and frequency of the sync
channel"):
  1. Sync times S_k = coarse + TDC fine (in coarse ticks).
  2. Index the sync records (detecting any missing records), and smooth the
     timing jitter with a moving linear fit (window `smooth_window` records).
     The pulse train is locked to the FPGA's 40 MHz clock, so its phase only
     drifts slowly relative to the MATRIX clock.
  3. Laser pulses between two sync records are equally spaced:
         t(32k + m) = S_k + m * (S_{k+1} - S_k) / 32.
  4. Each photon (after optional per-element delay correction) is assigned to
     the latest pulse at or before it: macrotime = pulse index, microtime =
     phase * laser period.

The time unit of the PTU is the laser period, which the metadata gives as
40 MHz (200 MHz base / 5, internal FPGA clock).  The MATRIX coarse tick length
in seconds is derived from it and reported.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class PulseTrain:
    sync_times: np.ndarray     # smoothed sync times, ticks (float64), one per record
    gap_pos: np.ndarray        # record positions k where sync records are missing before k
    gap_cum: np.ndarray        # cumulative number of missing records up to and incl. gap_pos
    divider: int               # laser pulses per sync record
    ticks_per_sync: float      # mean
    jitter_ticks_rms: float    # per-record jitter estimate before smoothing
    n_missing_sync: int

    def ticks_per_pulse(self):
        return self.ticks_per_sync / self.divider

    def sync_index(self, k):
        """Index (in sync periods) of record k, accounting for missing records."""
        if self.gap_pos.size == 0:
            return np.asarray(k, dtype=np.int64)
        j = np.searchsorted(self.gap_pos, k, side="right") - 1
        add = np.where(j >= 0, self.gap_cum[np.maximum(j, 0)], 0)
        return k + add

    def locate(self, t: np.ndarray):
        """Return (pulse_index:int64, phase:float in [0,1)) for times t (ticks).

        Times before the first / after the last sync record are extrapolated
        with the nearest period; the caller can count/flag them using
        `outside_mask`.
        """
        S = self.sync_times
        k = np.searchsorted(S, t, side="right") - 1
        k = np.clip(k, 0, S.size - 2)
        return self._locate(t, k, S[k], S[k + 1])

    def _locate(self, t, k, Sk, Sk1):
        """Pulse/phase of times t, given the sync record k at or before each t
        and the smoothed sync times Sk = S[k], Sk1 = S[k+1]."""
        i0 = self.sync_index(k)
        i1 = self.sync_index(k + 1)
        pulses = (i1 - i0) * self.divider          # divider unless a record is missing
        x = (t - Sk) / (Sk1 - Sk) * pulses
        m = np.floor(x)
        phase = x - m
        pulse = i0 * self.divider + m.astype(np.int64)
        return pulse, phase

    def outside_mask(self, t):
        return (t < self.sync_times[0]) | (t > self.sync_times[-1])

    def release(self):
        """Free the per-record sync times once nothing needs to be located."""
        self.sync_times = None


def _moving_mean(r: np.ndarray, half: int, block: int = 4_000_000) -> np.ndarray:
    """Centred moving mean (window 2*half+1) of r; at the edges the available
    one-sided window is used.  Computed in blocks to limit memory."""
    n = r.size
    w = 2 * half + 1
    out = np.empty_like(r)
    if n < w:
        out[:] = r.mean()
        return out
    for s0 in range(0, n, block):
        s1 = min(n, s0 + block)
        a = max(0, s0 - half)
        b = min(n, s1 + half)
        c = np.concatenate([[0.0], np.cumsum(r[a:b])])
        i = np.arange(s0, s1)
        lo = np.maximum(i - half, 0) - a
        hi = np.minimum(i + half + 1, n) - a
        out[s0:s1] = (c[hi] - c[lo]) / (hi - lo)
    return out


def build_pulse_train(sync_coarse, sync_fine, divider=32, smooth_window=64,
                      max_jitter_fraction=0.25) -> PulseTrain:
    S = sync_coarse.astype(np.float64)
    S += sync_fine
    return pulse_train_from_times(S, divider, smooth_window, max_jitter_fraction)


def pulse_train_from_times(S: np.ndarray, divider=32, smooth_window=64,
                           max_jitter_fraction=0.25, block: int = 4_000_000) -> PulseTrain:
    """Build the pulse train from sync times S (ticks, float64).

    S is modified in place and becomes PulseTrain.sync_times (no copies of the
    full array are made; all temporaries are block-sized).
    """
    n = S.size
    if n < 3:
        raise ValueError("fewer than 3 sync records; cannot build the pulse train")
    med = float(np.median(np.diff(S[: min(n, 1_000_001)])))
    dev, gap_pos, gap_cum, n_missing = _interval_scan(
        (S[s0:s0 + block] for s0 in range(0, n, block)), med, max_jitter_fraction)
    pt = PulseTrain(S, gap_pos, gap_cum, divider, 0.0, 0.0, n_missing)
    last = int(pt.sync_index(np.array([n - 1]))[0])
    # straight line through the ends; smooth the residual from it with a
    # centred moving mean, then replace S by line + smoothed residual
    S0 = float(S[0])
    slope = (S[-1] - S0) / last

    def resid(a, b):
        k = np.arange(a, b)
        return S[a:b] - (S0 + slope * pt.sync_index(k))

    kk = min(n, 2_000_000)
    jitter = float(np.diff(resid(0, kk), 2).std() / np.sqrt(6))
    half = max(1, smooth_window // 2)
    _smooth_inplace(S, resid, half, block)
    pt.jitter_ticks_rms = jitter
    pt.ticks_per_sync = float((S[-1] - S[0]) / last)
    return pt


def _interval_scan(chunks, med, max_jitter_fraction):
    """Check the sync intervals and find missing sync records.  `chunks` are
    consecutive pieces of the raw sync time array (any split gives the same
    result).  Returns (max deviation, gap_pos, gap_cum, n_missing)."""
    dev = 0.0
    gaps, gap_steps = [], []
    prev, off = None, 0          # last value of the previous chunk; index of chunk[0]
    for c in chunks:
        if prev is None:
            d, left0 = np.diff(c), off
        else:
            d, left0 = np.diff(c, prepend=prev), off - 1
        prev, off = c[-1], off + c.size
        if d.size == 0:
            continue
        steps = np.rint(d / med).astype(np.int32)
        if (steps < 1).any():
            raise ValueError("sync records closer together than expected - is the sync element correct?")
        d /= np.maximum(steps, 1)
        dev = max(dev, float(np.abs(d - med).max()))
        g = np.flatnonzero(steps > 1)
        gaps.append(g + left0)
        gap_steps.append(steps[g])
    if dev > max_jitter_fraction * med:
        raise ValueError(
            f"sync intervals are irregular (max deviation {dev:.2f} ticks of {med:.2f}); "
            "the sync is not a stable pulse train")
    gaps = np.concatenate(gaps) if gaps else np.empty(0, np.int64)
    gap_pos = (gaps + 1).astype(np.int64)                 # record index after the gap
    gs = np.concatenate(gap_steps) if gap_steps else np.empty(0, np.int32)
    gap_cum = np.cumsum(gs.astype(np.int64) - 1)
    n_missing = int(gap_cum[-1]) if gap_cum.size else 0
    return dev, gap_pos, gap_cum, n_missing


def _smoothed_block(raw, resid_line, n, half, block, j, i_end=None, chunk=1 << 19):
    """Smoothed sync times of records j*block .. i_end-1 (default: end of block j).

    Same arithmetic as _smooth_inplace() for that block: the moving-mean
    cumulative sum starts at the block's left halo.  It is computed in chunks;
    each chunk's cumulative sum starts from the previous chunk's last value, so
    the additions happen in the same order as in one np.cumsum (bit-identical).
    A partial block (i_end < end) also gives bit-identical values, because a
    prefix of a cumulative sum does not depend on later values."""
    s0 = j * block
    s1 = min(n, s0 + block)
    i_end = s1 if i_end is None else i_end
    a = max(0, s0 - half)
    b = min(n, i_end + half)
    # c[x] = sum of residuals a .. a+x-1  (c[0] = 0)
    c = np.empty(b - a + 1)
    c[0] = 0.0
    for u in range(a, b, chunk):
        v = min(b, u + chunk)
        r = raw(u, v) - resid_line(u, v)
        c[u - a + 1:v - a + 1] = np.cumsum(np.concatenate([[c[u - a]], r]))[1:]
    out = np.empty(i_end - s0)
    for u in range(s0, i_end, chunk):
        v = min(i_end, u + chunk)
        i = np.arange(u, v)
        lo = np.maximum(i - half, 0) - a
        hi = np.minimum(i + half + 1, n) - a
        sm = (c[hi] - c[lo]) / (hi - lo)
        S = raw(u, v)
        r = S - resid_line(u, v)
        o = out[u - s0:v - s0]
        np.subtract(S, r, out=o)
        o += sm
    return out


class StreamedPulseTrain(PulseTrain):
    """The same pulse train as pulse_train_from_times(), without holding the
    per-record sync times in memory: they are recomputed from the .raw file,
    one block of `block` records at a time, when photons or markers are
    located (a few blocks are cached).  Results are bit-identical to the
    in-memory version with the same block size."""

    CACHE_BLOCKS = 3

    def __init__(self, raw, n, divider, gap_pos, gap_cum, n_missing, half, block):
        super().__init__(None, gap_pos, gap_cum, divider, 0.0, 0.0, n_missing)
        self.n, self.half, self.block = n, half, block
        self._raw = raw
        self.n_blocks = -(-n // block)
        self._cache = {}

    def _setup_line(self, S0, slope):
        self._S0, self._slope = S0, slope

    def _line(self, a, b):
        return self._S0 + self._slope * self.sync_index(np.arange(a, b))

    def _smoothed(self, j, i_end=None):
        return _smoothed_block(self._raw, self._line, self.n, self.half, self.block, j, i_end)

    def _get(self, j):
        if j not in self._cache:
            if len(self._cache) >= self.CACHE_BLOCKS:
                self._cache.pop(next(iter(self._cache)))
            self._cache[j] = self._smoothed(j)
        return self._cache[j]

    def _value(self, i):
        """Smoothed time of single records i (rare edge cases only)."""
        return np.array([self._get(int(x) // self.block)[int(x) % self.block] for x in i])

    def locate(self, t: np.ndarray):
        """As PulseTrain.locate.  Times are grouped by block, and each block is
        computed once per call (blocks are visited in increasing order)."""
        t = np.asarray(t, dtype=np.float64)
        pulse = np.empty(t.size, np.int64)
        phase = np.empty(t.size, np.float64)
        if t.size == 0:
            return pulse, phase
        jb = np.clip(np.searchsorted(self.block_first, t, side="right") - 1, 0, self.n_blocks - 1)
        js = np.unique(jb)
        for j in js:
            m = slice(None) if js.size == 1 else (jb == j)
            tj = t[m]
            blk = self._get(int(j))
            s0 = int(j) * self.block
            k = np.clip(s0 + np.searchsorted(blk, tj, side="right") - 1, 0, self.n - 2)
            kl = k - s0
            inside = (kl >= 0) & (kl < blk.size)
            Sk = np.empty(tj.size)
            Sk1 = np.empty(tj.size)
            Sk[inside] = blk[kl[inside]]
            # S[k+1]: in this block, or the next block's first value at the block end
            k1 = kl[inside] + 1
            at_end = k1 >= blk.size
            nxt = self.block_first[j + 1] if j + 1 < self.n_blocks else np.nan
            Sk1[inside] = np.where(at_end, nxt, blk[np.minimum(k1, blk.size - 1)])
            if not inside.all():            # clipped into a neighbouring block
                o = ~inside
                Sk[o] = self._value(k[o])
                Sk1[o] = self._value(k[o] + 1)
            pulse[m], phase[m] = self._locate(tj, k, Sk, Sk1)
        return pulse, phase

    def outside_mask(self, t):
        return (t < self.first) | (t > self.last)

    def release(self):
        self._cache.clear()


def streamed_pulse_train(sync, table, divider=32, smooth_window=64,
                         max_jitter_fraction=0.25, block: int = 4_000_000) -> PulseTrain:
    """Pulse train from a rawstream.SyncReader without loading all sync records.

    `table` is the sync channel's TDC table (code -> fraction of a tick).
    Bit-identical to pulse_train_from_times(coarse + fine) with the same block.
    """
    n = sync.n
    half = max(1, smooth_window // 2)

    def raw(a, b):
        c, k = sync.read(a, b)
        S = c.astype(np.float64)
        f = table[k]
        if np.isnan(f).any():
            raise ValueError("sync element: TDC codes without calibration value")
        S += f
        return S

    # raw values around each block start, captured during the scan pass, so the
    # smoothed block-start times need no extra reading of the file
    n_blocks = -(-n // block)
    win = [(max(0, j * block - half), min(n, j * block + 1 + half)) for j in range(n_blocks)]
    win_val = [np.empty(b - a) for a, b in win]
    pos = [0]

    def capture(S):
        off = pos[0]
        for (a, b), v in zip(win, win_val):
            lo, hi = max(a, off), min(b, off + S.size)
            if lo < hi:
                v[lo - a:hi - a] = S[lo - off:hi - off]
        pos[0] = off + S.size

    def chunks():
        for c, k in sync.iter_blocks():
            S = c.astype(np.float64)
            f = table[k]
            if np.isnan(f).any():
                raise ValueError("sync element: TDC codes without calibration value")
            S += f
            capture(S)
            yield S

    if n < max(3, 2 * half + 1):        # tiny file: the in-memory version is simplest
        return pulse_train_from_times(raw(0, n), divider, smooth_window, max_jitter_fraction, block)

    med = float(np.median(np.diff(raw(0, min(n, 1_000_001)))))
    dev, gap_pos, gap_cum, n_missing = _interval_scan(chunks(), med, max_jitter_fraction)
    pt = StreamedPulseTrain(raw, n, divider, gap_pos, gap_cum, n_missing, half, block)
    last = int(pt.sync_index(np.array([n - 1]))[0])
    S0 = float(raw(0, 1)[0])
    slope = (raw(n - 1, n)[0] - S0) / last
    pt._setup_line(S0, slope)
    kk = min(n, 2_000_000)
    pt.jitter_ticks_rms = float(np.diff(raw(0, kk) - pt._line(0, kk), 2).std() / np.sqrt(6))
    # smoothed time of the first record of every block (needs only the halo)
    first = []
    for j, ((a, b), v) in enumerate(zip(win, win_val)):
        def raw_win(x, y, a=a, v=v):
            return v[x - a:y - a]
        first.append(_smoothed_block(raw_win, pt._line, n, half, block, j, j * block + 1)[0])
    pt.block_first = np.array(first)
    pt.first = pt.block_first[0]
    pt.last = pt._get(pt.n_blocks - 1)[-1]
    pt.ticks_per_sync = float((pt.last - pt.first) / last)
    sync.clear_cache()
    return pt


def _smooth_inplace(S, resid, half, block):
    """S <- S - r + moving_mean(r), r = resid(a, b) (the residual from the
    reference line).  Same arithmetic as _moving_mean() on the full residual
    array, but in place: the original values needed by the next block (its
    left halo) are saved before a block is overwritten."""
    n = S.size
    if n < 2 * half + 1:
        r = resid(0, n)
        S -= r
        S += r.mean()
        return
    saved = None            # original S values just before the current block
    for s0 in range(0, n, block):
        s1 = min(n, s0 + block)
        a = max(0, s0 - half)
        b = min(n, s1 + half)
        r = resid(a, b)
        if saved is not None:
            # S[a:s0] was overwritten by the previous block: use the saved originals
            r[: s0 - a] = saved
        c = np.concatenate([[0.0], np.cumsum(r)])
        i = np.arange(s0, s1)
        lo = np.maximum(i - half, 0) - a
        hi = np.minimum(i + half + 1, n) - a
        sm = (c[hi] - c[lo]) / (hi - lo)
        # residual of the next block's left halo, from original values
        nxt = max(s0, s1 - half)
        saved = r[nxt - a: s1 - a].copy()
        rb = r[s0 - a: s1 - a]
        S[s0:s1] -= rb
        S[s0:s1] += sm
