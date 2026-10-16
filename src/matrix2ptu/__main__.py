"""Command line interface.

    python -m matrix2ptu inspect  inputs/Matrix.raw [--msr inputs/Matrix.msr]
    python -m matrix2ptu convert  configs/matrix.yaml [--validate] [--tif inputs/Intensity.tif]
"""
from __future__ import annotations

import argparse
import json
import sys


def main(argv=None):
    ap = argparse.ArgumentParser(prog="matrix2ptu",
                                 description="Convert Abberior MATRIX .raw streams to PicoQuant PTU")
    sub = ap.add_subparsers(dest="cmd", required=True)

    a = sub.add_parser("inspect", help="summarise a .raw file (and its .msr metadata)")
    a.add_argument("raw")
    a.add_argument("--msr")

    c = sub.add_parser("convert", help="convert according to a YAML configuration")
    c.add_argument("config")
    c.add_argument("--validate", action="store_true", help="run validation checks afterwards")
    c.add_argument("--tif", help="reference intensity TIFF for validation (single-sequence data)")
    c.add_argument("--no-flopa", action="store_true", help="skip the tttrkit/FLOPA reconstruction check")
    c.add_argument("--dry-run", action="store_true",
                   help="decode and check parameters, write only the report (no PTU)")

    args = ap.parse_args(argv)
    try:
        if args.cmd == "inspect":
            from . import rawstream
            ev = rawstream.decode(args.raw)
            print(json.dumps(rawstream.summarize(ev), indent=2))
            if args.msr:
                from .msr_meta import read_msr_metadata
                m = read_msr_metadata(args.msr)
                print(json.dumps(m, indent=2, default=str))
        elif args.cmd == "convert":
            from .converter import convert
            res = convert(args.config, only_validate_inputs=args.dry_run)
            if args.validate and not args.dry_run:
                from .validate import validate
                out = validate(res, tif_path=args.tif, run_flopa=not args.no_flopa)
                if any(ch["passed"] is False for ch in out["checks"]):
                    return 2
    except Exception as e:  # concise error for users; full trace with MATRIX2PTU_DEBUG=1
        import os
        if os.environ.get("MATRIX2PTU_DEBUG"):
            raise
        print(f"ERROR: {e}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
