"""Minimal reader for Imspector .obf / .msr files (OMAS_BF container).

Used here only to read *metadata* (stack names/shapes and the Imspector XML
parameter documents) that help interpret the MATRIX .raw stream. It is not an
.obf-to-PTU converter.

Layout follows the publicly described OBF format (file header + linked list of
stacks, each with a fixed 368-byte header, optional zlib-compressed data and
a footer containing a tag dictionary). Only the parts needed here are parsed;
XML documents are located by scanning for ``<root ...>...</root>`` blocks.
"""
from __future__ import annotations

import re
import struct
import zlib
from dataclasses import dataclass, field

import numpy as np

FILE_MAGIC = b"OMAS_BF\n\xff\xff"
STACK_MAGIC = b"OMAS_BF_STACK\n\xff\xff"

# OBF data type codes -> numpy dtypes
OBF_DTYPES = {
    0x01: np.uint8, 0x02: np.int8, 0x04: np.uint16, 0x08: np.int16,
    0x10: np.uint32, 0x20: np.int32, 0x40: np.float32, 0x80: np.float64,
    0x400: np.uint64, 0x800: np.int64, 0x1000: np.bool_,
}


@dataclass
class ObfStack:
    index: int
    offset: int
    name: str
    description: str
    shape: tuple            # numpy order (slowest axis first)
    lengths_m: tuple        # physical length per axis (OBF order), metres or seconds
    offsets_m: tuple
    dtype: object
    compression: int
    data_offset: int
    data_len_disk: int
    xml_docs: list = field(default_factory=list)

    def read_data(self, buf: bytes) -> np.ndarray:
        raw = buf[self.data_offset:self.data_offset + self.data_len_disk]
        if self.compression:
            raw = zlib.decompress(raw)
        arr = np.frombuffer(raw, dtype=self.dtype)
        return arr.reshape(self.shape)


def _find_xml(blob: bytes) -> list[str]:
    docs = []
    for m in re.finditer(rb"<root[ >].*?</root>", blob, re.S):
        docs.append(m.group(0).decode("utf-8", "replace"))
    # UTF-16LE documents (used for some file-level descriptors)
    for m in re.finditer("<root[ >]".encode("utf-16-le") + rb".*?" +
                         "</root>".encode("utf-16-le"), blob, re.S):
        docs.append(m.group(0).decode("utf-16-le", "replace"))
    return docs


def read_obf(path: str):
    """Return (stacks, file_level_xml_docs, raw_bytes)."""
    with open(path, "rb") as f:
        buf = f.read()
    if not buf.startswith(FILE_MAGIC):
        raise ValueError(f"{path}: not an OBF/MSR file (bad magic)")
    (version,) = struct.unpack_from("<I", buf, 10)
    (first_stack,) = struct.unpack_from("<Q", buf, 14)

    stacks = []
    pos = first_stack
    idx = 0
    while pos:
        if buf[pos:pos + 16] != STACK_MAGIC:
            raise ValueError(f"{path}: bad stack magic at {pos}")
        p = pos + 16
        sver, rank = struct.unpack_from("<II", buf, p); p += 8
        res = struct.unpack_from("<15I", buf, p); p += 60
        lens = struct.unpack_from("<15d", buf, p); p += 120
        offs = struct.unpack_from("<15d", buf, p); p += 120
        dt, comp, _clevel, name_len, descr_len = struct.unpack_from("<5I", buf, p); p += 20
        _reserved, data_len, next_pos = struct.unpack_from("<QQQ", buf, p); p += 24
        name = buf[p:p + name_len].decode("utf-8", "replace"); p += name_len
        descr = buf[p:p + descr_len].decode("utf-8", "replace"); p += descr_len
        data_off = p
        end = next_pos if next_pos else len(buf)
        st = ObfStack(idx, pos, name, descr,
                      tuple(reversed(res[:rank])), lens[:rank], offs[:rank],
                      OBF_DTYPES.get(dt, np.uint8), comp, data_off, data_len)
        st.xml_docs = _find_xml(buf[data_off + data_len:end])
        stacks.append(st)
        pos = next_pos
        idx += 1
    file_docs = _find_xml(buf[:first_stack])
    return stacks, file_docs, buf
