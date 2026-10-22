"""Minimal reader for Imspector .obf / .msr files (OMAS_BF container).

Two kinds of OBF file occur with the MATRIX detector (see docs/OBF_REPORT.md):

* the measurement file (.msr, or an image .obf saved by Imspector/Lightbox):
  image / histogram stacks plus the Imspector parameter tree (XML);
* the *streamed* .obf written during acquisition when Matrix1 streaming is set
  to the OBF format (`/hwr/Matrix1/data/streaming/format` = 0; 4 = .raw): one
  1-D uint32 stack holding the same 32-bit event words as a .raw file, plus
  the parameter tree in that stack's footer.

Layout follows the publicly described OBF format (file header + linked list of
stacks, each with a fixed header, optional zlib-compressed data and a footer).
Only the parts needed here are parsed; XML documents are located by scanning
the footers for ``<root ...>...</root>`` blocks.  Headers and footers are read
with seeks, so large streamed files are never loaded whole.
"""
from __future__ import annotations

import re
import struct
import zlib
from dataclasses import dataclass, field

import numpy as np

FILE_MAGIC = b"OMAS_BF\n\xff\xff"
STACK_MAGIC = b"OMAS_BF_STACK\n\xff\xff"
STACK_FIXED_HEADER = 16 + 8 + 60 + 120 + 120 + 20 + 24   # bytes before name/description

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
    dtype_code: int
    compression: int
    data_offset: int
    data_len_disk: int
    path: str = ""
    xml_docs: list = field(default_factory=list)

    def read_data(self, buf: bytes | None = None) -> np.ndarray:
        """Stack contents as an array.  `buf` (the whole file) is optional;
        without it the data are read from the file."""
        if buf is not None:
            raw = buf[self.data_offset:self.data_offset + self.data_len_disk]
        else:
            with open(self.path, "rb") as f:
                f.seek(self.data_offset)
                raw = f.read(self.data_len_disk)
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


def iter_stacks(path: str, read_footers: bool = True):
    """Yield (ObfStack, footer_bytes) for every stack, reading only headers
    and (if read_footers) footers.  File-level XML documents: see read_obf()."""
    with open(path, "rb") as f:
        head = f.read(22)
        if not head.startswith(FILE_MAGIC):
            raise ValueError(f"{path}: not an OBF/MSR file (bad magic)")
        (first_stack,) = struct.unpack_from("<Q", head, 14)
        f.seek(0, 2)
        size = f.tell()
        pos, idx = first_stack, 0
        while pos:
            f.seek(pos)
            h = f.read(STACK_FIXED_HEADER)
            if h[:16] != STACK_MAGIC:
                raise ValueError(f"{path}: bad stack magic at {pos}")
            p = 16
            _sver, rank = struct.unpack_from("<II", h, p); p += 8
            res = struct.unpack_from("<15I", h, p); p += 60
            lens = struct.unpack_from("<15d", h, p); p += 120
            offs = struct.unpack_from("<15d", h, p); p += 120
            dt, comp, _clevel, name_len, descr_len = struct.unpack_from("<5I", h, p); p += 20
            _reserved, data_len, next_pos = struct.unpack_from("<QQQ", h, p); p += 24
            nd = f.read(name_len + descr_len)
            name = nd[:name_len].decode("utf-8", "replace")
            descr = nd[name_len:].decode("utf-8", "replace")
            data_off = pos + STACK_FIXED_HEADER + name_len + descr_len
            st = ObfStack(idx, pos, name, descr,
                          tuple(reversed(res[:rank])), lens[:rank], offs[:rank],
                          OBF_DTYPES.get(dt, np.uint8), dt, comp, data_off, data_len, path)
            footer = b""
            if read_footers:
                end = next_pos if next_pos else size
                f.seek(data_off + data_len)
                footer = f.read(max(0, end - data_off - data_len))
                st.xml_docs = _find_xml(footer)
            yield st, footer
            pos = next_pos
            idx += 1


def read_obf(path: str):
    """Return (stacks, file_level_xml_docs, None).  The third element used to
    be the whole file; ObfStack.read_data() now reads from the file itself."""
    stacks = [st for st, _ in iter_stacks(path)]
    with open(path, "rb") as f:
        head = f.read(22)
        (first_stack,) = struct.unpack_from("<Q", head, 14)
        f.seek(0)
        file_docs = _find_xml(f.read(first_stack))
    return stacks, file_docs, None


@dataclass
class EventStream:
    """Location of the MATRIX event words inside a streamed .obf file."""
    path: str
    stack_name: str
    byte_offset: int
    n_words: int            # including any trailing zero padding


def find_event_stream(path: str) -> EventStream:
    """Locate the stack holding the raw MATRIX event words in a streamed .obf.

    Accepted: exactly one stack that is 1-D, uint32 and uncompressed.  Anything
    else (image stacks only, several candidate stacks, compressed data) is
    refused with an explanation rather than guessed."""
    stacks = [st for st, _ in iter_stacks(path, read_footers=False)]
    cand = [s for s in stacks if len(s.shape) == 1 and s.dtype_code == 0x10]
    if not cand:
        desc = ", ".join(f"'{s.name}' {s.shape} {np.dtype(s.dtype).name}" for s in stacks)
        raise ValueError(
            f"{path}: no 1-D uint32 event-stream stack found (stacks: {desc}). This looks like an "
            "image/histogram OBF (e.g. a saved measurement), not a streamed MATRIX acquisition; "
            "it holds no photon events and cannot be converted to event-based PTU.")
    if len(cand) > 1:
        raise ValueError(f"{path}: {len(cand)} candidate event-stream stacks "
                         f"({[s.name for s in cand]}); not supported (not seen in test data)")
    s = cand[0]
    if s.compression:
        raise ValueError(f"{path}: event-stream stack '{s.name}' is compressed; not supported "
                         "(not seen in test data)")
    if s.data_len_disk % 4 or s.data_len_disk != 4 * s.shape[0]:
        raise ValueError(f"{path}: stack '{s.name}' length {s.data_len_disk} bytes does not match "
                         f"{s.shape[0]} uint32 words")
    return EventStream(path, s.name, s.data_offset, s.shape[0])
