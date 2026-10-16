"""Writer for PicoQuant .ptu files with event-based T3 records.

Record format: PicoQuant "HydraHarp V2 / MultiHarp T3" 32-bit record
(TTResultFormat_TTTRRecType = 0x01010304, the same code that Imspector writes
for MultiHarp data, see the SPAD.ptu reference):

    bit 31      special
    bits 30..25 channel (6 bit)   photons: detector channel 0..62
                                  special: 63 = overflow, 1..15 = marker bits
    bits 24..10 dtime   (15 bit)  micro time in units of MeasDesc_Resolution
    bits  9..0  nsync   (10 bit)  macro time in sync (laser) periods, mod 1024
                                  overflow: number of 1024-wraps represented

Header: "PQTTTR\\0\\0" + version + tag list (PicoQuant unified tag format:
32-byte identifier, int32 index, uint32 type code, 8-byte value / length,
followed by data for strings) terminated by Header_End.

This module only *encodes*.  Correctness is checked by reading the files back
with independent readers (ptufile, tttrkit) in validate.py.
"""
from __future__ import annotations

import datetime as _dt
import struct
import uuid

import numpy as np

REC_TYPE_HYDRAHARP2_T3 = 0x01010304
WRAP = 1024
OVERFLOW_CHANNEL = 63

# tag type codes
TY_EMPTY8 = 0xFFFF0008
TY_BOOL8 = 0x00000008
TY_INT8 = 0x10000008
TY_FLOAT8 = 0x20000008
TY_TDATETIME = 0x21000008
TY_ANSISTRING = 0x4001FFFF


def _tag(ident: str, value, index: int = -1, typ: int | None = None) -> bytes:
    if len(ident) > 31:
        raise ValueError(f"tag identifier too long: {ident}")
    head = ident.encode("ascii").ljust(32, b"\0")
    if typ is None:
        if value is None:
            typ = TY_EMPTY8
        elif isinstance(value, bool):
            typ = TY_BOOL8
        elif isinstance(value, (int, np.integer)):
            typ = TY_INT8
        elif isinstance(value, (float, np.floating)):
            typ = TY_FLOAT8
        elif isinstance(value, _dt.datetime):
            typ = TY_TDATETIME
        elif isinstance(value, str):
            typ = TY_ANSISTRING
        else:
            raise TypeError(f"unsupported tag value type for {ident}: {type(value)}")
    if typ == TY_EMPTY8:
        return head + struct.pack("<iIq", index, typ, 0)
    if typ == TY_BOOL8:
        return head + struct.pack("<iIq", index, typ, 1 if value else 0)
    if typ == TY_INT8:
        return head + struct.pack("<iIq", index, typ, int(value))
    if typ == TY_FLOAT8:
        return head + struct.pack("<iId", index, typ, float(value))
    if typ == TY_TDATETIME:
        days = (value - _dt.datetime(1899, 12, 30)).total_seconds() / 86400.0
        return head + struct.pack("<iId", index, typ, days)
    if typ == TY_ANSISTRING:
        data = value.encode("latin-1", "replace") + b"\0"
        data = data.ljust((len(data) + 7) // 8 * 8, b"\0")
        return head + struct.pack("<iIq", index, typ, len(data)) + data
    raise ValueError(typ)


def build_header(tags: list[tuple]) -> tuple[bytes, int]:
    """tags: list of (ident, value[, index]).  Returns (bytes, offset of the
    NumberOfRecords value so it can be patched later)."""
    parts = [b"PQTTTR\0\0", b"1.0.00\0\0"]
    pos = 16
    nrec_off = -1
    for t in tags:
        ident, value = t[0], t[1]
        index = t[2] if len(t) > 2 else -1
        b = _tag(ident, value, index)
        if ident == "TTResult_NumberOfRecords":
            nrec_off = pos + 40
        parts.append(b)
        pos += len(b)
    parts.append(_tag("Header_End", None))
    if nrec_off < 0:
        raise ValueError("header must contain TTResult_NumberOfRecords")
    return b"".join(parts), nrec_off


def encode_t3(nsync: np.ndarray, channel: np.ndarray, dtime: np.ndarray,
              special: np.ndarray, hi_prev: int = 0) -> np.ndarray:
    """Encode time-sorted events (absolute nsync, int64) into T3 records,
    inserting overflow records.  Returns uint32 array.  `hi_prev` is the
    number of wraps already emitted (for block-wise encoding)."""
    nsync = np.asarray(nsync, dtype=np.int64)
    if nsync.size and (np.diff(nsync) < 0).any():
        raise ValueError("events are not sorted by nsync")
    if nsync.size and nsync[0] < 0:
        raise ValueError("negative nsync")
    if (np.asarray(dtime) >= 1 << 15).any() or (np.asarray(dtime) < 0).any():
        raise ValueError("dtime out of 15-bit range")
    if (np.asarray(channel) > 63).any():
        raise ValueError("channel out of 6-bit range")

    hi = nsync // WRAP
    lo = (nsync % WRAP).astype(np.uint32)
    # number of wraps to emit before each event
    prev = np.concatenate([[0], hi[:-1]])
    need = hi - prev
    if need.size:
        need[0] = hi[0] - hi_prev
    # Each overflow record carries up to 1023 wraps; events needing more get
    # several overflow records.
    n_ovf = np.where(need > 0, (need + 1022) // 1023, 0)
    total = nsync.size + int(n_ovf.sum())
    out = np.empty(total, dtype=np.uint32)
    ev_pos = np.arange(nsync.size) + np.cumsum(n_ovf)
    ev_words = ((special.astype(np.uint32) & 1) << 31) | \
               ((channel.astype(np.uint32) & 0x3F) << 25) | \
               ((dtime.astype(np.uint32) & 0x7FFF) << 10) | lo
    out[ev_pos] = ev_words
    # fill overflow records
    idx = np.flatnonzero(n_ovf)
    for i in idx:  # usually few (one per ~25 us of silence); loop is fine
        n = int(need[i])
        start = ev_pos[i] - n_ovf[i]
        k = 0
        while n > 0:
            c = min(n, 1023)
            out[start + k] = (1 << 31) | (OVERFLOW_CHANNEL << 25) | c
            n -= c
            k += 1
    return out


def encode_t3_fast(nsync, channel, dtime, special, hi_prev: int = 0):
    """Vectorised variant of encode_t3 for the common case where no gap
    exceeds 1023 wraps (always true here: markers/photons are dense, and
    sync-derived nsync gaps of >26 ms would be needed). Falls back otherwise."""
    nsync = np.asarray(nsync, dtype=np.int64)
    hi = nsync // WRAP
    need = np.diff(hi, prepend=np.int64(hi_prev))
    if (need > 1023).any() or (need < 0).any():
        return encode_t3(nsync, channel, dtime, special, hi_prev)
    has = need > 0
    total = nsync.size + int(has.sum())
    out = np.empty(total, dtype=np.uint32)
    ev_pos = np.arange(nsync.size) + np.cumsum(has)
    out[ev_pos] = ((special.astype(np.uint32) & 1) << 31) | \
                  ((channel.astype(np.uint32) & 0x3F) << 25) | \
                  ((dtime.astype(np.uint32) & 0x7FFF) << 10) | \
                  (nsync % WRAP).astype(np.uint32)
    ovf_pos = ev_pos[has] - 1
    out[ovf_pos] = (np.uint32(1) << 31) | (np.uint32(OVERFLOW_CHANNEL) << 25) | \
                   need[has].astype(np.uint32)
    return out


def write_ptu(path: str, tags: list[tuple], records: np.ndarray) -> None:
    header, nrec_off = build_header(tags)
    with open(path, "wb") as f:
        f.write(header)
        records.astype("<u4", copy=False).tofile(f)
        f.seek(nrec_off)
        f.write(struct.pack("<q", int(records.size)))


def write_ptu_t3(path: str, tags: list[tuple], nsync, channel, dtime, special,
                 block: int = 4_000_000) -> int:
    """Encode time-sorted events and write them, block by block (the full
    record array is never held in memory).  Byte-identical to
    write_ptu(path, tags, encode_t3_fast(...)).  Returns the record count."""
    header, nrec_off = build_header(tags)
    nrec = 0
    hi_prev = 0
    with open(path, "wb") as f:
        f.write(header)
        for s0 in range(0, len(nsync), block):
            sl = slice(s0, s0 + block)
            rec = encode_t3_fast(nsync[sl], channel[sl], dtime[sl], special[sl], hi_prev)
            rec.astype("<u4", copy=False).tofile(f)
            nrec += rec.size
            hi_prev = int(nsync[sl][-1]) // WRAP
        f.seek(nrec_off)
        f.write(struct.pack("<q", nrec))
    return nrec


def new_guid() -> str:
    return "{" + str(uuid.uuid4()).upper() + "}"
