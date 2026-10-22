"""Command line interface.

    python -m matrix2ptu inspect  inputs/Matrix.raw [--msr inputs/Matrix.msr]
    python -m matrix2ptu inspect  inputs2/MatrixOBF.obf [--msr inputs2/MatrixOBF.msr]
    python -m matrix2ptu convert  configs/matrix.yaml [--validate] [--tif inputs/Intensity.tif]
    python -m matrix2ptu ism-calibrate --seq CONF=beads_seq0.ptu --seq STED=beads_seq1.ptu \
        --detect-in STED --metadata Beads.obf --out calibration.json
    python -m matrix2ptu ism-reassign target_seq0.ptu --calibration calibration.json --sequence CONF \
        --out-dir outdir [--metadata Target.obf]
"""
from __future__ import annotations

import argparse
import json
import sys


def main(argv=None):
    ap = argparse.ArgumentParser(prog="matrix2ptu",
                                 description="Convert Abberior MATRIX .raw / streamed .obf event streams to PicoQuant PTU")
    sub = ap.add_subparsers(dest="cmd", required=True)

    a = sub.add_parser("inspect", help="summarise a .raw or streamed .obf file (and its metadata)")
    a.add_argument("raw", help=".raw stream or streamed .obf")
    a.add_argument("--msr")

    c = sub.add_parser("convert", help="convert according to a YAML configuration")
    c.add_argument("config")
    c.add_argument("--validate", action="store_true", help="run validation checks afterwards")
    c.add_argument("--tif", help="reference intensity TIFF for validation (single-sequence data)")
    c.add_argument("--no-flopa", action="store_true", help="skip the tttrkit/FLOPA reconstruction check")
    c.add_argument("--dry-run", action="store_true",
                   help="decode and check parameters, write only the report (no PTU)")

    k = sub.add_parser("ism-calibrate", help="bead calibration of per-element pixel-reassignment shifts")
    k.add_argument("--seq", action="append", required=True, metavar="NAME=PTU",
                   help="sequence name and its per-sequence detector-resolved PTU (repeat per sequence)")
    k.add_argument("--detect-in", required=True, help="sequence whose pooled image is used to find beads (STED)")
    k.add_argument("--metadata", help=".obf/.msr of the bead acquisition (stored as fingerprint)")
    k.add_argument("--out", required=True, help="calibration JSON")
    k.add_argument("--isolation-nm", type=float, default=320)

    r = sub.add_parser("ism-reassign", help="apply a shift calibration to a per-sequence detector-resolved PTU")
    r.add_argument("ptu")
    r.add_argument("--calibration", required=True)
    r.add_argument("--sequence", required=True, help="calibration sequence to use (e.g. CONF)")
    r.add_argument("--out-dir", required=True)
    r.add_argument("--metadata", help=".obf/.msr of the target (settings compared with the calibration)")
    r.add_argument("--y-rounding", choices=("nearest", "stochastic"), default="nearest")
    r.add_argument("--no-verify", action="store_true")

    args = ap.parse_args(argv)
    try:
        if args.cmd == "inspect":
            from . import rawstream
            ev = rawstream.decode(args.raw)
            print(json.dumps(rawstream.summarize(ev), indent=2))
            if rawstream.is_obf(args.raw):
                from .msr_meta import read_msr_metadata
                print("embedded .obf parameters:")
                print(json.dumps(read_msr_metadata(args.raw), indent=2, default=str))
            if args.msr:
                from .msr_meta import read_msr_metadata
                m = read_msr_metadata(args.msr)
                print(json.dumps(m, indent=2, default=str))
        elif args.cmd == "convert":
            from .converter import convert
            res = convert(args.config, only_validate_inputs=args.dry_run,
                          keep_validation_data=args.validate)
            if args.validate and not args.dry_run:
                from .validate import validate
                out = validate(res, tif_path=args.tif, run_flopa=not args.no_flopa)
                if any(ch["passed"] is False for ch in out["checks"]):
                    return 2
        elif args.cmd == "ism-calibrate":
            from . import ism
            seqs = dict(a.split("=", 1) for a in args.seq)
            cal = ism.calibrate(seqs, args.detect_in, args.metadata, args.out, isolation_nm=args.isolation_nm)
            for n, c in cal["sequences"].items():
                print(f"{n}: {len(c['beads_used'])} beads; {c['recommendation']}")
            print(f"wrote {args.out}")
        elif args.cmd == "ism-reassign":
            import os
            from . import ism
            cal = ism.load_calibration(args.calibration)
            os.makedirs(args.out_dir, exist_ok=True)
            stem = os.path.splitext(os.path.basename(args.ptu))[0]
            stem = stem.replace("_original", "").replace("detector-resolved", "").replace("detector_resolved", "")
            stem = stem.strip("_-")
            det = os.path.join(args.out_dir, f"{stem}_detector-resolved_reassigned.ptu")
            pool = os.path.join(args.out_dir, f"{stem}_pooled_reassigned.ptu")
            rep = ism.reassign_ptu(args.ptu, cal, args.sequence, det, pool, calibration_file=args.calibration,
                                   target_metadata=args.metadata, y_rounding=args.y_rounding)
            for w in rep["warnings"]:
                print("WARNING:", w)
            if not args.no_verify:
                rep["checks"] = ism.verify_reassigned(rep)
            rep.pop("_arrays")
            with open(os.path.join(args.out_dir, f"{stem}_reassignment_report.json"), "w") as f:
                json.dump(ism._jsonable(rep), f, indent=1)
            print(json.dumps(ism._jsonable(rep["stats"]), indent=1))
            print(f"wrote {det}, {pool}")
    except Exception as e:  # concise error for users; full trace with MATRIX2PTU_DEBUG=1
        import os
        if os.environ.get("MATRIX2PTU_DEBUG"):
            raise
        print(f"ERROR: {e}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
