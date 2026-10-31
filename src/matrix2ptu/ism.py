"""Pixel reassignment (ISM) workflow for MATRIX data converted by matrix2ptu.

    calibrate():  bead calibration of per-element shifts, one table per
                  acquisition sequence, saved as JSON (shifts in nm)
    reassign_ptu(): apply a calibration to one per-sequence detector-resolved PTU,
                  writing reassigned detector-resolved and pooled PTUs
    verify_reassigned(): structural / content checks of a reassigned PTU

Inputs are the per-sequence detector-resolved PTU files written with
``output.split_sequences: true`` (one PTU per acquisition sequence; the line
accumulation of the sequence is read from the UsrMX_Sequences tag).

Shift convention (also written into every calibration file):
  x = scan (fast) axis, + to the right = later in the line; y = line index,
  + downwards; the same orientation as Imspector's image stacks.
  shift s_e of element e = position of the element image relative to the
  pooled (all-element) image: I_e(r) ~ I_ref(r - s_e).  Reassignment moves every
  photon of element e by -s_e, so the reassigned image stays registered to the
  original pooled image.
"""
from __future__ import annotations

import datetime as dt
import json
import os
import re

import numpy as np

from . import beadcal
from .reassign import (N_EL, load_photons, merge_events, reassign_positions, reassign_times,
                       write_derived_ptu)
from .validate import iter_ptu_records

CAL_FORMAT = "matrix2ptu ISM shift calibration v1"


# --------------------------------------------------------------------------- PTU scan info
def scan_info(ptu_path: str) -> dict:
    import ptufile
    p = ptufile.PtuFile(ptu_path)
    try:
        t = dict(p.tags)
    finally:
        p.close()
    seq = str(t.get("UsrMX_Sequences", "0:seq:x1"))
    parts = seq.split(";")
    if len(parts) != 1:
        raise ValueError(f"{ptu_path}: contains {len(parts)} sequences ({seq}); reassignment needs one PTU "
                         "per sequence (convert with output.split_sequences: true)")
    m = re.match(r"\d+:(.*):x(\d+)$", parts[0])
    if not m:
        raise ValueError(f"{ptu_path}: cannot parse UsrMX_Sequences '{seq}'")
    name, acc = m.group(1), int(m.group(2))
    mode = str(t.get("UsrMX_OutputMode", ""))
    if mode != "detector_resolved":
        raise ValueError(f"{ptu_path}: output mode '{mode}'; reassignment needs the detector-resolved PTU "
                         "(the pooled PTU has lost the element identity of each photon)")
    return {"sequence": name, "accumulation": acc, "pixels_x": int(t["ImgHdr_PixX"]),
            "pixels_y": int(t["ImgHdr_PixY"]) // acc, "pixel_size_nm": float(t["ImgHdr_PixResol"]) * 1e3,
            "tags": t}


def open_sequence(ptu_path: str):
    si = scan_info(ptu_path)
    pd = load_photons(ptu_path, si["accumulation"], si["pixels_x"], si["pixels_y"], si["pixel_size_nm"])
    return si, pd


# --------------------------------------------------------------------------- calibration
SEQ_DEFAULTS = {  # fit half-window and start sigma (pixels at 20 nm); scaled by pixel size
    "confocal": {"half_win_nm": 240, "sigma_nm": 120},
    "sted": {"half_win_nm": 120, "sigma_nm": 42},
}


def calibrate(seq_ptus: dict[str, str], sted_sequence: str, metadata_source: str | None,
              out_json: str, kinds: dict[str, str] | None = None, isolation_nm: float = 320,
              max_neighbours: int = 1) -> dict:
    """Bead calibration.  seq_ptus = {sequence name: detector-resolved PTU};
    beads are detected in the pooled image of `sted_sequence` (best separation).
    kinds = {sequence name: 'confocal' | 'sted'} (default: name contains STED)."""
    from .msr_meta import optical_fingerprint
    images, infos = {}, {}
    for name, path in seq_ptus.items():
        si, pd = open_sequence(path)
        images[name] = pd.images()
        infos[name] = {k: v for k, v in si.items() if k != "tags"} | {"file": os.path.basename(path),
                                                                         "photons_in_image": int(pd.inside.sum())}
    det_img = images[sted_sequence].sum(0)
    cal = {"format": CAL_FORMAT, "created": dt.datetime.now().isoformat(timespec="seconds"),
           "conventions": {
               "axes": "x = scan (fast) axis, + right (later within the line); y = line index, + down; "
                       "orientation of Imspector's image stacks",
               "shift": "s_e = position of element e's image relative to the pooled image: I_e(r) ~ I_ref(r - s_e)",
               "reassignment": "move each photon of element e by -s_e",
               "elements": "data element 0..22 = PTU channel 0..22 = Imspector layout label 1..23",
               "units": "shift_nm; divide by the target pixel size for pixels"},
           "fingerprint": optical_fingerprint(metadata_source) if metadata_source else None,
           "sequences": {}}
    for name, stack in images.items():
        kind = (kinds or {}).get(name, "sted" if "STED" in name.upper() else "confocal")
        px = infos[name]["pixel_size_nm"]
        d = SEQ_DEFAULTS[kind]
        res = beadcal.calibrate_sequence(stack, det_img, px, int(round(d["half_win_nm"] / px)),
                                         d["sigma_nm"] / px, isolation_px=isolation_nm / px,
                                         max_neighbours=max_neighbours)
        loo = res["loo_validation"]["summary"]["k1"]["fwhm_ratio"]
        beneficial = loo["median"] < 0.99 and loo["n_improved"] >= 0.8 * loo["n"]
        res["recommendation"] = ("apply: held-out beads narrow consistently" if beneficial else
                                 "no measurable benefit on held-out beads (reassignment optional)")
        res["kind"] = kind
        res["input"] = infos[name]
        res["shift_nm"] = res["shift_px"] * px
        res["uncertainty_nm"] = res["uncertainty_px"] * px
        cal["sequences"][name] = res
    save_calibration(cal, out_json)
    return cal


def _jsonable(v):
    if isinstance(v, np.ndarray):
        return v.tolist()
    if isinstance(v, (np.floating, np.integer, np.bool_)):
        return v.item()
    if isinstance(v, dict):
        return {k: _jsonable(x) for k, x in v.items()}
    if isinstance(v, (list, tuple)):
        return [_jsonable(x) for x in v]
    return v


def save_calibration(cal: dict, path: str):
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(_jsonable(cal), f, indent=1)


def load_calibration(path: str) -> dict:
    with open(path, encoding="utf-8") as f:
        cal = json.load(f)
    if cal.get("format") != CAL_FORMAT:
        raise ValueError(f"{path}: not a {CAL_FORMAT} file")
    return cal


# --------------------------------------------------------------------------- reassignment
def reassign_ptu(src_ptu: str, cal: dict, sequence: str, out_detres: str, out_pooled: str | None,
                 calibration_file: str = "", target_metadata: str | None = None, note: str = "",
                 y_rounding: str = "nearest", seed: int = 0) -> dict:
    """Write reassigned PTU(s) for one sequence.  Returns a report dict and keeps
    the arrays needed for verification under the key '_arrays'."""
    from .msr_meta import compare_fingerprints, optical_fingerprint
    si, pd = open_sequence(src_ptu)
    if sequence not in cal["sequences"]:
        raise ValueError(f"calibration has no sequence '{sequence}' (has {list(cal['sequences'])})")
    cs = cal["sequences"][sequence]
    shifts_nm = np.asarray(cs["shift_nm"], float)
    shifts_px = shifts_nm / si["pixel_size_nm"]
    warnings = []
    if si["sequence"] != sequence:
        warnings.append(f"PTU sequence '{si['sequence']}' reassigned with calibration of '{sequence}'")
    fp_diff = None
    if target_metadata and cal.get("fingerprint"):
        fp_diff = compare_fingerprints(cal["fingerprint"], optical_fingerprint(target_metadata))
        for d in fp_diff:
            if d["relevance"] not in ("none", "info"):
                warnings.append(f"setting differs from the calibration acquisition: {d['key']} "
                                f"({d['calibration']} -> {d['target']}; relevance: {d['relevance']})")
    tn, keep, stats, dy = reassign_times(pd, shifts_px, y_rounding=y_rounding, seed=seed)
    tags = [
        ("UsrMX_Reassigned", True),
        ("UsrMX_ReassignCalibration", os.path.basename(calibration_file)[:200]),
        ("UsrMX_ReassignSequence", sequence),
        ("UsrMX_ReassignMethod", f"per-element shift -s_e: x exact to 1 laser period, y to whole lines "
                                 f"({y_rounding} rounding); photons leaving the field dropped; micro time and "
                                 "channel unchanged"),
        ("UsrMX_ReassignMacroTime", "macro time of moved photons encodes the new position, not the detection time"),
        ("UsrMX_SourcePTU", os.path.basename(src_ptu)),
    ]
    for e in range(N_EL):
        tags.append(("UsrMX_ReassignShiftX_nm", float(shifts_nm[e, 0]), e))
        tags.append(("UsrMX_ReassignShiftY_nm", float(shifts_nm[e, 1]), e))
    comment = (f"Pixel-reassigned ({sequence}) from {os.path.basename(src_ptu)} with "
               f"{os.path.basename(calibration_file)}. Moved photons carry a position-coded macro time; "
               f"micro time and detector channel are the original values. {note}").strip()
    files = []
    for path, pooled in ((out_detres, False), (out_pooled, True)):
        if not path:
            continue
        n, c, d, s = merge_events(pd, tn, keep, pooled)
        t2 = tags + [("UsrMX_OutputMode", "pooled_reassigned" if pooled else "detector_resolved_reassigned")]
        nrec = write_derived_ptu(src_ptu, path, t2, n, c, d, s, replace_comment=comment)
        files.append({"file": path, "pooled": pooled, "records": int(nrec), "photons": int(keep.sum())})
    return {"source": src_ptu, "sequence": sequence, "scan": {k: v for k, v in si.items() if k != "tags"},
            "shifts_px": shifts_px, "stats": stats, "files": files, "warnings": warnings,
            "fingerprint_differences": fp_diff,
            "_arrays": {"pd": pd, "tn": tn, "keep": keep, "shifts_px": shifts_px, "dy": dy}}


# --------------------------------------------------------------------------- verification
def _ptu_image(path, pixels_x, pixels_y, acc, channels):
    ls, le = [], []
    for _, d, t in iter_ptu_records(path):
        ls.append(t[(d["marker"] & 1) > 0]); le.append(t[(d["marker"] & 2) > 0])
    ls, le = np.concatenate(ls), np.concatenate(le)
    img = np.zeros((channels, pixels_y, pixels_x), np.int64)
    dt_by_ch = [[] for _ in range(channels)]
    n_ph, n_mk, mono = 0, 0, True
    last = -1
    for _, d, t in iter_ptu_records(path):
        ev = (d["channel"] >= 0) | (d["marker"] > 0)
        te = t[ev]
        if te.size and (te[0] < last or (np.diff(te) < 0).any()):
            mono = False
        if te.size:
            last = te[-1]
        ph = (d["channel"] >= 0) & (d["marker"] == 0)
        n_mk += int((d["marker"] > 0).sum())
        tp, ch, dtm = t[ph], d["channel"][ph].astype(np.int64), d["dtime"][ph]
        n_ph += tp.size
        for c in range(channels):
            dt_by_ch[c].append(dtm[ch == c])
        k = np.searchsorted(ls, tp, side="right") - 1
        kk = np.clip(k, 0, None)
        ok = (k >= 0) & (tp < le[kk])
        x = np.floor((tp - ls[kk]) / (le[kk] - ls[kk]) * pixels_x).astype(np.int64)
        y = kk // acc
        ok &= (x >= 0) & (x < pixels_x) & (y < pixels_y)
        np.add.at(img, (ch[ok], y[ok], x[ok]), 1)
    return img, [np.concatenate(v) for v in dt_by_ch], {"photons": n_ph, "markers": n_mk, "monotonic": mono,
                                                         "line_starts": int(ls.size), "line_stops": int(le.size)}


def verify_reassigned(rep: dict, run_tttrkit: bool = True) -> list[dict]:
    """Checks for each reassigned file of one reassign_ptu() result."""
    a = rep["_arrays"]
    pd, tn, keep, shifts, dy = a["pd"], a["tn"], a["keep"], a["shifts_px"], a["dy"]
    si = rep["scan"]
    nx, ny, acc = si["pixels_x"], si["pixels_y"], si["accumulation"]
    # expected image from the photon table (independent of the PTU encoding)
    e = pd.ch.astype(np.int64)
    m = pd.inside & keep
    xe = np.floor(pd.x[m] - shifts[e[m], 0]).astype(np.int64)
    ye = pd.y[m] - dy[m]
    exp_img = np.zeros((N_EL, ny, nx), np.int64)
    np.add.at(exp_img, (e[m], ye, xe), 1)
    # same, but from the new macro times (tests the PTU encoding exactly)
    tk = tn[keep]
    k = np.searchsorted(pd.ls, tk, side="right") - 1
    kk = np.clip(k, 0, None)
    ok = (k >= 0) & (tk < pd.le[kk])
    xt = np.floor((tk - pd.ls[kk]) / (pd.le[kk] - pd.ls[kk]) * nx).astype(np.int64)
    yt = kk // acc
    ok &= (xt >= 0) & (xt < nx) & (yt < ny)
    exp_t = np.zeros((N_EL, ny, nx), np.int64)
    np.add.at(exp_t, (e[keep][ok], yt[ok], xt[ok]), 1)
    checks = []
    for f in rep["files"]:
        ch = 1 if f["pooled"] else N_EL
        img, dts, info = _ptu_image(f["file"], nx, ny, acc, ch)
        expect = exp_img.sum(0, keepdims=True) if f["pooled"] else exp_img
        expect_t = exp_t.sum(0, keepdims=True) if f["pooled"] else exp_t
        # vs floor(x - sx): the integer macro time can put a photon 1 px lower when
        # x - sx lies within one laser period (1/4000 px) above a pixel edge
        diff = int(np.abs(img - expect).sum())
        diff_t = int(np.abs(img - expect_t).sum())
        orig_dt = [np.sort(pd.dtime[keep & (pd.ch == c)]) for c in range(N_EL)]
        if f["pooled"]:
            orig_dt = [np.sort(pd.dtime[keep])]
        dt_same = all(np.array_equal(np.sort(a_), b_) for a_, b_ in zip(dts, orig_dt))
        per_ch = [int(v.size) for v in dts]
        exp_ch = [int(v.size) for v in orig_dt]
        c = {"file": os.path.basename(f["file"]),
             "structure": {"readable_ptufile": True, "time_monotonic": info["monotonic"],
                           "line_start_markers": info["line_starts"], "line_stop_markers": info["line_stops"],
                           "markers_same_as_source": info["line_starts"] == pd.ls.size and info["line_stops"] == pd.le.size},
             "photon_counts": {"file": info["photons"], "expected": int(keep.sum()),
                               "per_channel_equal": per_ch == exp_ch},
             "micro_time_multiset_per_channel_identical": bool(dt_same),
             "image_vs_new_macro_times": {"sum_abs_diff": diff_t, "identical": diff_t == 0},
             "image_vs_continuous_positions": {"sum_abs_diff": diff, "fraction": diff / max(int(expect.sum()), 1)}}
        if run_tttrkit:
            from .validate import flopa_reconstruct
            r = {"pixels_y": ny, "pixels_x": nx, "frames": 1,
                 "sequences": [{"line_accumulation": acc}]}
            ds, _ = flopa_reconstruct(f["file"], r, max_detector=64)
            pc = ds["photon_count"].values[0, 0]          # (line, pixel, channel)
            tk = np.moveaxis(pc, -1, 0)[:ch]
            c["tttrkit_reconstruction"] = {"shape": list(ds["photon_count"].shape),
                                           "identical_to_own": bool(np.array_equal(tk, img)),
                                           "sum_abs_diff": int(np.abs(tk - img).sum())}
        checks.append(c)
    return checks
