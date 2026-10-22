"""Control test: does reading the event stream out of a streamed .obf change anything?

The event words of the .obf are copied into a plain headerless .raw file in a
temporary directory.  Both are then converted with the same configuration
(the given OBF config, with input.obf replaced by input.raw for the copy, and
the .obf's own embedded parameters supplied via the .msr setting being kept).
The PTU record sections must be byte-identical and the headers may differ only
in provenance tags.

Usage:
    .venv/bin/python tools/obf_container_check.py configs/matrix_obf.yaml
"""
import copy
import os
import sys
import tempfile

import numpy as np
import yaml

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from matrix2ptu.converter import convert          # noqa: E402
from matrix2ptu.obf import find_event_stream      # noqa: E402

ALLOWED_TAG_DIFFS = {"File_GUID", "File_CreatingTime", "File_Comment", "UsrMX_SourceFile",
                     "UsrMX_SourceFormat", "UsrMX_LaserFreqSource"}


def records(path):
    b = open(path, "rb").read()
    i = b.index(b"Header_End") + 48          # tag: 32 ident + 4 idx + 4 type + 8 value
    return b[i:]


def tags(path):
    import ptufile
    with ptufile.PtuFile(path) as p:
        return dict(p.tags)


def main(cfg_path):
    cfg = yaml.safe_load(open(cfg_path))
    obf = cfg["input"]["obf"]
    es = find_event_stream(obf)
    with tempfile.TemporaryDirectory() as tmp:
        raw = os.path.join(tmp, os.path.splitext(os.path.basename(obf))[0] + ".raw")
        words = np.fromfile(obf, "<u4", count=es.n_words, offset=es.byte_offset)
        words.tofile(raw)
        print(f"copied {es.n_words} words from stack '{es.stack_name}' (offset {es.byte_offset}) to {raw}")
        results = {}
        for kind in ("obf", "raw"):
            c = copy.deepcopy(cfg)
            c["output"]["directory"] = os.path.join(tmp, kind)
            c["output"]["split_sequences"] = False
            if kind == "raw":
                c["input"] = {"raw": raw, "msr": cfg["input"].get("msr")}
                # the .raw route needs metadata from a .msr; to compare like with like the
                # embedded .obf parameters must equal the .msr (checked by the obf run)
            p = os.path.join(tmp, f"{kind}.yaml")
            yaml.safe_dump(c, open(p, "w"))
            results[kind] = convert(p)
        rep = results["obf"]["report"]
        if rep["parameters"].get("msr_vs_obf_parameters", {}).get("value") != "agree":
            print("NOTE: embedded .obf parameters and .msr were not confirmed equal; "
                  "differences below may come from the metadata, not the container")
        ok = True
        for fo, fr in zip(results["obf"]["report"]["outputs"], results["raw"]["report"]["outputs"]):
            same_rec = records(fo["file"]) == records(fr["file"])
            to, tr = tags(fo["file"]), tags(fr["file"])
            diff = sorted(k for k in set(to) | set(tr) if to.get(k) != tr.get(k))
            bad = [k for k in diff if k not in ALLOWED_TAG_DIFFS]
            print(f"{fo['mode']}: records byte-identical = {same_rec} ({fo['records']} records); "
                  f"header tags differing: {diff}; unexpected: {bad}")
            ok &= same_rec and not bad
        print("CONTAINER CHECK", "PASS" if ok else "FAIL")
        return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1]))
