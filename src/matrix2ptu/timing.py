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
        i0 = self.sync_index(k)
        i1 = self.sync_index(k + 1)
        pulses = (i1 - i0) * self.divider          # divider unless a record is missing
        x = (t - S[k]) / (S[k + 1] - S[k]) * pulses
        m = np.floor(x)
        phase = x - m
        pulse = i0 * self.divider + m.astype(np.int64)
        return pulse, phase

    def outside_mask(self, t):
        return (t < self.sync_times[0]) | (t > self.sync_times[-1])


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


def coarse_to_float_inplace(a: np.ndarray, block: int = 4_000_000) -> np.ndarray:
    """Convert an int64 array to float64 in its own memory (exact below 2**53;
    unwrapped coarse times are far below).  The int64 array must not be used
    afterwards.  Returns the float64 view."""
    if a.dtype != np.int64:
        raise TypeError("expected int64")
    f = a.view(np.float64)
    for s0 in range(0, a.size, block):
        f[s0:s0 + block] = a[s0:s0 + block].astype(np.float64)
    return f


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
    # interval check and gap detection, block-wise
    dev = 0.0
    gaps, gap_steps = [], []
    for s0 in range(0, n - 1, block):
        s1 = min(n - 1, s0 + block)
        d = np.diff(S[s0:s1 + 1])
        steps = np.rint(d / med).astype(np.int32)
        if (steps < 1).any():
            raise ValueError("sync records closer together than expected - is the sync element correct?")
        d /= np.maximum(steps, 1)
        dev = max(dev, float(np.abs(d - med).max()))
        g = np.flatnonzero(steps > 1)
        gaps.append(g + s0)
        gap_steps.append(steps[g])
    if dev > max_jitter_fraction * med:
        raise ValueError(
            f"sync intervals are irregular (max deviation {dev:.2f} ticks of {med:.2f}); "
            "the sync is not a stable pulse train")
    gaps = np.concatenate(gaps)
    gap_pos = (gaps + 1).astype(np.int64)                 # record index after the gap
    gap_cum = np.cumsum(np.concatenate(gap_steps).astype(np.int64) - 1)
    n_missing = int(gap_cum[-1]) if gap_cum.size else 0
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
