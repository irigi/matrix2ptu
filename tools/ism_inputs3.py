"""Pixel-reassignment round on inputs3 (Beads = calibration, Dense = target).

    .venv/bin/python tools/ism_inputs3.py all          # every stage below, in order
    .venv/bin/python tools/ism_inputs3.py convert      # .obf -> per-sequence PTU (+ validated combined PTU)
    .venv/bin/python tools/ism_inputs3.py calibrate    # bead calibration (Beads)
    .venv/bin/python tools/ism_inputs3.py transfer     # does the Beads calibration fit Dense?
    .venv/bin/python tools/ism_inputs3.py reassign     # reassigned PTUs + checks, both datasets
    .venv/bin/python tools/ism_inputs3.py validate     # image / FLIM / FLOPA comparisons
    .venv/bin/python tools/ism_inputs3.py figures

Everything is written below outputs/inputs3_reassignment/ (figures to docs/figures/).
"""
import glob
import json
import os
import shutil
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "src"))
from matrix2ptu import beadcal as bc               # noqa: E402
from matrix2ptu import ism                         # noqa: E402
from matrix2ptu.msr_meta import read_flim_stacks   # noqa: E402

R = os.path.join(ROOT, "outputs", "inputs3_reassignment")
INP = os.path.join(ROOT, "inputs3")
DATASETS = ("Beads", "Dense")
SEQS = (("CONF", 1), ("STED", 3))
CAL_FULL = os.path.join(R, "2_calibration", "Beads_ism_calibration.json")
CAL_COMPACT = os.path.join(ROOT, "configs", "calibration", "MATRIX1_ism_shifts_2026-10-22_Beads640_8AU.json")
PY = os.path.join(ROOT, ".venv", "bin", "python")


def name(D, i, s, mode, state):
    return f"{D}_seq{i}-{s}_{mode}_{state}.ptu"


def orig(D, i, s, mode):
    return os.path.join(R, "1_original", D, name(D, i, s, mode, "original"))


def reas(D, i, s, mode):
    return os.path.join(R, "3_reassigned", D, name(D, i, s, mode, "reassigned"))


def dump(obj, path):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        json.dump(ism._jsonable(obj), f, indent=1)
    print("wrote", os.path.relpath(path, ROOT))


# =========================================================================== convert
def stage_convert():
    for D in DATASETS:
        lc = D.lower()
        for cfg, extra in ((f"configs/{lc}_obf.yaml", ["--validate"]), (f"configs/{lc}_obf_split.yaml", [])):
            print("converting", cfg)
            subprocess.run([PY, "-m", "matrix2ptu", "convert", cfg, *extra], cwd=ROOT, check=True,
                           env=dict(os.environ, PYTHONPATH=os.path.join(ROOT, "src")))
        d = os.path.join(R, "1_original", D)
        for i, (s, _) in enumerate(SEQS):
            for mode, m2 in (("detector_resolved", "detector-resolved"), ("pooled", "pooled")):
                src = os.path.join(d, f"{D}_{mode}_seq{i}_{s}.ptu")
                if os.path.exists(src):
                    os.replace(src, orig(D, i, s, m2))


# =========================================================================== calibrate
def compact(cal):
    out = {k: cal[k] for k in ("format", "created", "conventions", "fingerprint")}
    out["sequences"] = {}
    for n, c in cal["sequences"].items():
        out["sequences"][n] = {
            "kind": c["kind"], "recommendation": c["recommendation"],
            "pixel_size_nm_at_calibration": c["pixel_size_nm"],
            "shift_nm": c["shift_nm"], "uncertainty_nm": c["uncertainty_nm"],
            "n_beads": len(c["beads_used"]),
            "loo_fwhm_ratio_k1_median": c["loo_validation"]["summary"]["k1"]["fwhm_ratio"]["median"],
            "source": c["input"]["file"]}
    return out


def stage_calibrate():
    ptus = {s: orig("Beads", i, s, "detector-resolved") for i, (s, _) in enumerate(SEQS)}
    cal = ism.calibrate(ptus, "STED", os.path.join(INP, "Beads.obf"), CAL_FULL)
    ism.save_calibration(compact(cal), CAL_COMPACT)
    print("wrote", os.path.relpath(CAL_FULL, ROOT), "and", os.path.relpath(CAL_COMPACT, ROOT))
    for n, c in cal["sequences"].items():
        print(n, "beads", len(c["beads_used"]), "max |shift|",
              round(float(np.hypot(*np.asarray(c["shift_nm"]).T).max()), 1), "nm;", c["recommendation"])


# =========================================================================== shared helpers
_cache = {}


def element_images(D, i, s):
    key = (D, i, s)
    if key not in _cache:
        si, pd = ism.open_sequence(orig(D, i, s, "detector-resolved"))
        rng = np.random.default_rng(12345)
        half = rng.random(pd.t.size) < 0.5          # binomial split -> independent Poisson halves
        _cache[key] = (pd.images(), pd.images(sel=half), pd.images(sel=~half))
    return _cache[key]


def apr_shifts(S):
    """Self-estimated shifts: cross-correlation of each element with the sum of the
    OTHER elements (avoids the zero-shift noise self-correlation), converted to the
    pooled-image reference (s_pooled = s_LOO * (1 - photon fraction))."""
    tot = S.sum(0)
    w = S.reshape(len(S), -1).sum(1) / S.sum()
    return np.array([bc.xcorr_shift(S[e], tot - S[e], win=12) for e in range(len(S))]) * (1 - w)[:, None]


def pattern_scale(ref, S):
    A, _, res = bc.fit_affine(ref, S, np.ones(len(ref)))
    d = bc.describe_affine(A)
    return {"scale": d["scale_per_pitch"], "rotation_deg": d["rotation_deg_x_right_y_down"],
            "residual_rms_px": float(np.sqrt((res ** 2).sum(1).mean()))}


BANDS = [(0.02, 0.05), (0.05, 0.08), (0.08, 0.11), (0.11, 0.14), (0.14, 0.18)]


def band_gain(h0, h1, shifts_list, crop=10):
    """Amplitude gain per spatial-frequency band of the pooled image after
    reassignment with each shift table, relative to no reassignment.  Signal power
    per band = cross-power of two independent photon halves (noise-free in
    expectation); significance from the FRC of the unreassigned halves."""
    H = h0.shape[1] - 2 * crop
    win = np.hanning(H)[:, None] * np.hanning(H)[None, :]
    yy, xx = np.indices((H, H))
    fr = np.hypot(yy - H // 2, xx - H // 2) / H

    def spec(img):
        img = img[crop:-crop, crop:-crop]
        return np.fft.fftshift(np.fft.fft2((img - img.mean()) * win))

    def pooled(S, sh):
        return sum(bc.fourier_shift(S[e], -sh[e, 0], -sh[e, 1]) for e in range(len(S)))

    A0, B0 = spec(h0.sum(0)), spec(h1.sum(0))
    rows = []
    for lo, hi in BANDS:
        m = (fr >= lo) & (fr < hi)
        c = (A0 * np.conj(B0)).real[m]
        frc = c.sum() / np.sqrt((np.abs(A0[m]) ** 2).sum() * (np.abs(B0[m]) ** 2).sum())
        rows.append({"band_period_nm": [20 / hi, 20 / lo], "frc_before": float(frc),
                     "significant": bool(frc > 3 / np.sqrt(m.sum() / 2)), "power_before": float(c.sum())})
    out = {"bands": rows, "gains": {}}
    for lab, sh in shifts_list.items():
        A, B = spec(pooled(h0, sh)), spec(pooled(h1, sh))
        g = []
        for (lo, hi), r in zip(BANDS, rows):
            m = (fr >= lo) & (fr < hi)
            p = (A * np.conj(B)).real[m].sum()
            g.append(float(np.sqrt(max(p, 0) / r["power_before"])) if r["power_before"] > 0 else None)
        out["gains"][lab] = g
    return out


def object_widths(D, s, hw, sg, R_iso=10):
    """Gaussian sigma (px) of quasi-isolated objects in the pooled image of D/s
    (neighbours fitted jointly), objects detected in D's STED image."""
    sted = element_images(D, 1, "STED")[0].sum(0)
    allsrc, _ = bc.detect_beads(sted, rel_thresh=0.0, noise_k=5)
    pk, _ = bc.detect_beads(sted, rel_thresh=0.05)
    cl = bc.classify_beads(pk, sted.shape, R_iso, hw + 2, sted, all_sources=allsrc, max_ellipticity=1.0)
    img = element_images(D, 0 if s == "CONF" else 1, s)[0].sum(0)
    sig = []
    for b in cl:
        if not b["isolated"]:
            continue
        y0, x0 = int(b["y"]), int(b["x"])
        w = img[y0 - hw:y0 + hw + 1, x0 - hw:x0 + hw + 1]
        nb = [(yy - y0 + hw, xx - x0 + hw) for yy, xx, _ in allsrc
              if 0.5 < np.hypot(yy - b["y"], xx - b["x"]) and max(abs(yy - y0), abs(xx - x0)) <= hw + 2 * sg]
        p, e, ok = bc.fit_gauss(w, b["y"] - y0 + hw, b["x"] - x0 + hw, sg, nb)
        if ok:
            sig.append(p[3])
    sig = np.array(sig)
    return {"n": int(sig.size), "sigma_px_median": float(np.median(sig)),
            "sigma_px_iqr": np.percentile(sig, [25, 75]).tolist(), "fwhm_nm_median": float(np.median(sig) * 2.3548 * 20)}


# =========================================================================== transfer
def stage_transfer():
    from matrix2ptu.msr_meta import compare_fingerprints, optical_fingerprint
    cal = ism.load_calibration(CAL_FULL)
    out = {"fingerprint_differences": compare_fingerprints(cal["fingerprint"],
                                                           optical_fingerprint(os.path.join(INP, "Dense.obf")))}
    for i, (s, _) in enumerate(SEQS):
        c = np.asarray(cal["sequences"][s]["shift_px"])
        r = {}
        for D in DATASETS:
            full, h0, h1 = element_images(D, i, s)
            a = apr_shifts(full)
            q = []
            for ys in (slice(0, 125), slice(125, 250)):
                for xs in (slice(0, 125), slice(125, 250)):
                    q.append(pattern_scale(c, apr_shifts(full[:, ys, xs]))["scale"])
            r[D] = {"self_estimate_vs_bead_calibration": pattern_scale(c, a),
                    "self_estimate_quadrant_scales": q,
                    "self_estimate_shift_px": a,
                    "band_gain": band_gain(h0, h1, {"bead calibration": c, "0.5 x bead": 0.5 * c,
                                                    "1.4 x bead": 1.4 * c, "self-estimate": a})}
        hw, sg = (12, 6.0) if s == "CONF" else (6, 2.1)
        for D in DATASETS:
            r[D]["object_widths"] = object_widths(D, s, hw, sg)
        out[s] = r
        print(s, {D: (round(r[D]["self_estimate_vs_bead_calibration"]["scale"], 3),
                      [None if g is None else round(g, 3) for g in r[D]["band_gain"]["gains"]["bead calibration"]],
                      round(r[D]["object_widths"]["fwhm_nm_median"], 1)) for D in DATASETS})
    # blur test: are the estimators object-size dependent? (Beads images convolved with a common kernel)
    from scipy import ndimage
    blur = {}
    for i, (s, _) in enumerate(SEQS):
        c = np.asarray(cal["sequences"][s]["shift_px"])
        full = element_images("Beads", i, s)[0]
        blur[s] = {}
        for sk in (0, 1, 2, 3, 4, 6):
            S = np.array([ndimage.gaussian_filter(full[e], sk) for e in range(23)]) if sk else full
            blur[s][f"sigma{sk}px"] = pattern_scale(c, apr_shifts(S))
    out["blur_test_self_estimate"] = blur
    # band-limited self-estimates: which shift aligns coarse vs fine detail?
    # (reliable only while the shift is well below half the band period: STED yes, CONF fine bands no)
    bl = {}
    for i, (s, _) in enumerate(SEQS):
        c = np.asarray(cal["sequences"][s]["shift_px"])
        bl[s] = {}
        for D in DATASETS:
            S = element_images(D, i, s)[0]
            w = S.reshape(23, -1).sum(1) / S.sum()
            ky = np.fft.fftfreq(S.shape[1])[:, None]
            kx = np.fft.fftfreq(S.shape[2])[None, :]
            f = np.hypot(ky, kx)
            bl[s][D] = {}
            for lo, hi in ((0.0, 0.05), (0.05, 0.08), (0.08, 0.11), (0.11, 0.18)):
                B = np.real(np.fft.ifft2(np.fft.fft2(S) * ((f >= lo) & (f < hi))))
                tot = B.sum(0)
                a = np.array([bc.xcorr_shift(B[e], tot - B[e], win=12) for e in range(23)]) * (1 - w)[:, None]
                bl[s][D][f"{20 / hi:.0f}-{20 / lo:.0f}nm" if lo else f">{20 / hi:.0f}nm"] = pattern_scale(c, a)
    out["band_limited_self_estimate"] = bl
    dump(out, os.path.join(R, "4_validation", "transfer_assessment.json"))


# =========================================================================== reassign
def stage_reassign():
    import tifffile
    cal = ism.load_calibration(CAL_FULL)
    for D in DATASETS:
        os.makedirs(os.path.join(R, "3_reassigned", D, "images"), exist_ok=True)
        reports = []
        for i, (s, acc) in enumerate(SEQS):
            note = ""
            if D == "Dense":
                note = "Shifts transferred from the Beads acquisition (see transfer assessment)."
            rep = ism.reassign_ptu(orig(D, i, s, "detector-resolved"), cal, s,
                                   reas(D, i, s, "detector-resolved"), reas(D, i, s, "pooled"),
                                   calibration_file=CAL_FULL, target_metadata=os.path.join(INP, f"{D}.obf"),
                                   note=note)
            checks = ism.verify_reassigned(rep)
            a = rep["_arrays"]
            pd = a["pd"]
            before = pd.images()
            e = pd.ch.astype(np.int64)
            m = pd.inside & a["keep"]
            after = np.zeros_like(before)
            np.add.at(after, (e[m], (pd.y[m] - a["dy"][m]).astype(np.int64),
                              np.floor(pd.x[m] - a["shifts_px"][e[m], 0]).astype(np.int64)), 1)
            sh = a["shifts_px"]
            exact = np.array([bc.fourier_shift(before[k], -sh[k, 0], -sh[k, 1]) for k in range(23)])
            base = os.path.join(R, "3_reassigned", D, "images", f"{D}_seq{i}-{s}")
            tifffile.imwrite(base + "_pooled_original.tif", before.sum(0).astype(np.float32))
            tifffile.imwrite(base + "_pooled_reassigned-events.tif", after.sum(0).astype(np.float32))
            tifffile.imwrite(base + "_pooled_reassigned-fourier.tif", exact.sum(0).astype(np.float32))
            tifffile.imwrite(base + "_elements_reassigned-events.tif", after.astype(np.float32),
                             metadata={"axes": "CYX"})
            rep.pop("_arrays")
            rep["checks"] = checks
            reports.append(rep)
            print(D, s, {k: rep["stats"][k] for k in ("photons_in_image", "photons_moved",
                                                      "photons_dropped_out_of_field")}, rep["warnings"])
            for c in checks:
                print("   ", c["file"], "monotonic", c["structure"]["time_monotonic"], "counts ok",
                      c["photon_counts"]["per_channel_equal"], "dtime identical",
                      c["micro_time_multiset_per_channel_identical"], "image identical",
                      c["image_vs_new_macro_times"]["identical"], "tttrkit identical",
                      c.get("tttrkit_reconstruction", {}).get("identical_to_own"))
        dump(reports, os.path.join(R, "3_reassigned", D, f"{D}_reassignment_report.json"))


# =========================================================================== detector / FRC diagnostics
def detector_correlations(D, i=0, s="CONF"):
    """Same-laser-pulse photon pairs between element pairs vs the accidental rate
    expected from the per-pixel counts (crosstalk shows up between neighbours)."""
    si, pd = ism.open_sequence(orig(D, i, s, "detector-resolved"))
    m = pd.inside
    t, e = pd.t[m], pd.ch[m].astype(np.int64)
    px = pd.y[m] * 250 + np.floor(pd.x[m]).astype(np.int64)
    o = np.argsort(t, kind="stable")
    t, e, px = t[o], e[o], px[o]
    same = np.flatnonzero(np.diff(t) == 0)
    C = np.zeros((23, 23))
    np.add.at(C, (e[same], e[same + 1]), 1)
    C = C + C.T
    imgs = np.zeros((23, 62500))
    np.add.at(imgs, (e, px), 1)
    pulses_per_pixel = float(np.median(pd.le - pd.ls)) / 250
    Exp = imgs @ imgs.T / pulses_per_pixel
    np.fill_diagonal(Exp, 0)
    P = bc.layout()
    dist = np.hypot(*(P[:, None, :] - P[None, :, :]).transpose(2, 0, 1))
    nb, far = (dist > 0.9) & (dist < 1.1), dist > 1.9
    return {"sequence": s, "photons": int(t.size),
            "neighbour_pairs_observed": float(C[nb].sum() / 2), "neighbour_pairs_accidental": float(Exp[nb].sum() / 2),
            "distant_pairs_observed": float(C[far].sum() / 2), "distant_pairs_accidental": float(Exp[far].sum() / 2),
            "excess_neighbour_pairs_fraction_of_photons": float((C[nb].sum() - Exp[nb].sum()) / 2 / t.size)}


FRC_BANDS = [(0.02, 0.05), (0.05, 0.08), (0.08, 0.11), (0.11, 0.14), (0.14, 0.18), (0.18, 0.24)]


def frc_multi(D, i, s, shifts, n_split=8):
    si, pd = ism.open_sequence(orig(D, i, s, "detector-resolved"))
    out = {"bands_period_nm": [[20 / hi, 20 / lo] for lo, hi in FRC_BANDS]}
    for k in (0.0, 1.0):
        F, X = [], []
        for seed in range(n_split):
            sel = np.random.default_rng(100 + seed).random(pd.t.size) < 0.5
            h0, h1 = pd.images(sel=sel), pd.images(sel=~sel)
            a = sum(bc.fourier_shift(h0[e], -k * shifts[e, 0], -k * shifts[e, 1]) for e in range(23))[10:240, 10:240]
            b = sum(bc.fourier_shift(h1[e], -k * shifts[e, 0], -k * shifts[e, 1]) for e in range(23))[10:240, 10:240]
            f, c = bc.frc(a, b)
            F.append([float(np.mean(c[(f >= lo) & (f < hi)])) for lo, hi in FRC_BANDS])
            X.append(float(bc.frc_resolution(f, c) * 20))
        out[f"k{k:g}"] = {"frc_mean": np.mean(F, 0).tolist(), "frc_sd": np.std(F, 0).tolist(),
                          "crossing_1_7_nm_mean": float(np.mean(X)), "crossing_1_7_nm_sd": float(np.std(X))}
    return out


# =========================================================================== validate
def ptu_image(path, acc, channels):
    from matrix2ptu.ism import _ptu_image
    return _ptu_image(path, 250, 250, acc, channels)


def bead_fwhm_table(D, i, s, images, beads, allsrc, hw, sg):
    rows = []
    for b in beads:
        y0, x0 = int(round(b["y"])), int(round(b["x"]))
        nb = [(yy - y0 + hw, xx - x0 + hw) for yy, xx, _ in allsrc
              if 0.5 < np.hypot(yy - b["y"], xx - b["x"]) and max(abs(yy - y0), abs(xx - x0)) <= hw + 2 * sg]
        r = {"x": b["x"], "y": b["y"]}
        for k, im in images.items():
            w = im[y0 - hw:y0 + hw + 1, x0 - hw:x0 + hw + 1]
            if w.shape != (2 * hw + 1, 2 * hw + 1):
                break
            p, e, ok = bc.fit_gauss(w, b["y"] - y0 + hw, b["x"] - x0 + hw, sg, nb)
            r[k] = {"sigma_px": float(p[3]), "amp": float(p[0]), "ok": ok}
        else:
            rows.append(r)
    return rows


def flopa_maps(path, acc):
    from napari_flopa.core.io.loader import read_ptu_file
    from napari_flopa.core.processing.reconstruction import reconstruct_ptu_to_dataset
    from tttrkit.ptuio.reconstructor import ScanConfig
    data = read_ptu_file(path, header=False)
    cfg = ScanConfig(lines=250, pixels=250, frames=1, line_accumulations=(acc,))
    ds = reconstruct_ptu_to_dataset(data, cfg, outputs=["photon_count", "mean_arrival_time", "phasor"])
    g = lambda k: ds[k].values[0, 0, :, :, 0].astype(float)
    return {"dims": dict(ds.sizes), "count": g("photon_count"), "mat": g("mean_arrival_time"),
            "g": g("phasor_g"), "s": g("phasor_s"), "res_ps": float(ds["tcspc_resolution"]) * 1e12}


def stage_validate():
    cal = ism.load_calibration(CAL_FULL)
    out = {}
    for D in DATASETS:
        cubes = read_flim_stacks(os.path.join(INP, f"{D}.msr"))
        cube_names = list(cubes)
        sted = element_images(D, 1, "STED")[0].sum(0)
        allsrc, _ = bc.detect_beads(sted, rel_thresh=0.0, noise_k=5)
        pk, _ = bc.detect_beads(sted, rel_thresh=0.05)
        for i, (s, acc) in enumerate(SEQS):
            key = f"{D}_{s}"
            r = {}
            # --- originals: detector-resolved vs pooled vs Imspector cube
            det, _, info_d = ptu_image(orig(D, i, s, "detector-resolved"), acc, 23)
            pool, _, info_p = ptu_image(orig(D, i, s, "pooled"), acc, 1)
            ref = cubes[cube_names[i]].sum(-1)
            r["original"] = {
                "detector_resolved_sum_equals_pooled": bool(np.array_equal(det.sum(0), pool[0])),
                "photons_in_image": int(pool.sum()),
                "line_markers_detres_pooled": [info_d["line_starts"], info_p["line_starts"]],
                "vs_imspector_cube": {"cube": cube_names[i], "counts_ptu": int(pool.sum()), "counts_imspector": int(ref.sum()),
                                      "identical_pixels_fraction": float((pool[0] == ref).mean()),
                                      "sum_abs_diff_fraction": float(np.abs(pool[0] - ref).sum() / ref.sum()),
                                      "r": float(np.corrcoef(pool[0].ravel(), ref.ravel())[0, 1])}}
            # --- after: event-level vs exact sub-pixel (image domain)
            aft, _, _ = ptu_image(reas(D, i, s, "pooled"), acc, 1)
            fou = __import__("tifffile").imread(os.path.join(R, "3_reassigned", D, "images",
                                                             f"{D}_seq{i}-{s}_pooled_reassigned-fourier.tif"))
            b = 8
            inner = (slice(b, -b), slice(b, -b))
            r["event_vs_fourier"] = {"r_inner": float(np.corrcoef(aft[0][inner].ravel(), fou[inner].ravel())[0, 1]),
                                     "photons_after": int(aft.sum())}
            # --- residual misalignment after reassignment (self-estimate on reassigned elements)
            aft_el, _, _ = ptu_image(reas(D, i, s, "detector-resolved"), acc, 23)
            sh_cal = np.asarray(cal["sequences"][s]["shift_px"])
            res_after = apr_shifts(aft_el[:, inner[0], inner[1]].astype(float))
            res_before = apr_shifts(det[:, inner[0], inner[1]].astype(float))
            r["self_estimate_shift_rms_px"] = {"before": float(np.sqrt((res_before ** 2).sum(1).mean())),
                                               "after": float(np.sqrt((res_after ** 2).sum(1).mean()))}
            # --- bead / object widths before vs after (event-level and exact)
            hw, sg = (12, 6.0) if s == "CONF" else (6, 2.1)
            imgs = {"before": pool[0].astype(float), "after_events": aft[0].astype(float), "after_fourier": fou}
            c_cal = cal["sequences"][s]
            if D == "Beads":
                sets = {"calibration_beads": c_cal["beads_used"],
                        "crowded_or_overlapping": [b_ for b_ in c_cal["bead_table"] if not b_["isolated"] and not b_["edge"]]}
            else:
                cl = bc.classify_beads(pk, sted.shape, 10, hw + 2, sted, all_sources=allsrc, max_ellipticity=1.0)
                sets = {"quasi_isolated_objects": [b_ for b_ in cl if b_["isolated"]]}
            r["object_fwhm"] = {}
            for nm, beads in sets.items():
                rows = bead_fwhm_table(D, i, s, imgs, beads, allsrc, hw, sg)
                if not rows:
                    continue
                rat = np.array([x["after_events"]["sigma_px"] / x["before"]["sigma_px"] for x in rows])
                pk_r = np.array([x["after_events"]["amp"] / x["before"]["amp"] for x in rows])
                r["object_fwhm"][nm] = {
                    "n": len(rows), "fwhm_before_nm_median": float(np.median([x["before"]["sigma_px"] for x in rows]) * 47.096),
                    "fwhm_after_nm_median": float(np.median([x["after_events"]["sigma_px"] for x in rows]) * 47.096),
                    "ratio_after_over_before_median": float(np.median(rat)),
                    "ratio_iqr": np.percentile(rat, [25, 75]).tolist(),
                    "n_narrower": int((rat < 1).sum()), "peak_ratio_median": float(np.median(pk_r)),
                    "ratio_fourier_median": float(np.median([x["after_fourier"]["sigma_px"] / x["before"]["sigma_px"] for x in rows])),
                    "per_object": rows}
            # --- edge effects: photons lost per pixel row/column near the border
            lost = (pool[0].astype(float) - aft[0].astype(float))
            r["edge"] = {"max_abs_shift_px": float(np.abs(sh_cal).max()),
                         "border_px_affected": int(np.ceil(np.abs(sh_cal).max())),
                         "photons_lost_fraction": float(lost.sum() / pool.sum())}
            # --- FLIM: pooled decay and FLOPA maps before vs after
            fl_b = flopa_maps(orig(D, i, s, "pooled"), acc)
            fl_a = flopa_maps(reas(D, i, s, "pooled"), acc)
            ok = (fl_b["count"] >= 100) & (fl_a["count"] >= 100)
            ok[:b] = ok[-b:] = False
            ok[:, :b] = ok[:, -b:] = False
            wsum = lambda fl, k: float(np.nansum(fl[k][ok] * fl["count"][ok]) / fl["count"][ok].sum())
            ns = fl_b["res_ps"] * 1e-3            # FLOPA/tttrkit mean arrival time is in raw TCSPC bins
            r["flopa"] = {"dims": fl_a["dims"], "images_equal_own_reconstruction":
                          bool(np.array_equal(fl_a["count"], aft[0].astype(float))),
                          "pixels_compared": int(ok.sum()),
                          "count_weighted_mean_arrival_ns": [wsum(fl_b, "mat") * ns, wsum(fl_a, "mat") * ns],
                          "count_weighted_phasor_g": [wsum(fl_b, "g"), wsum(fl_a, "g")],
                          "count_weighted_phasor_s": [wsum(fl_b, "s"), wsum(fl_a, "s")],
                          "per_pixel_mean_arrival_change_ns_median_abs":
                              float(np.nanmedian(np.abs(fl_a["mat"][ok] - fl_b["mat"][ok]))) * ns}
            out[key] = r
            print(key, json.dumps({k: v for k, v in r.items() if k not in ("object_fwhm",)}, default=str)[:900])
            for nm, v in r["object_fwhm"].items():
                print("   ", nm, {k: v[k] for k in v if k != "per_object"})
    out["detector_crosstalk"] = {D: detector_correlations(D) for D in DATASETS}
    out["frc_multi_split"] = {f"{D}_{s}": frc_multi(D, i, s, np.asarray(cal["sequences"][s]["shift_px"]))
                              for D in DATASETS for i, (s, _) in enumerate(SEQS)}
    print(json.dumps(out["detector_crosstalk"], indent=1))
    dump(out, os.path.join(R, "4_validation", "validation.json"))


if __name__ == "__main__":
    st = sys.argv[1] if len(sys.argv) > 1 else "all"
    stages = {"convert": stage_convert, "calibrate": stage_calibrate, "transfer": stage_transfer,
              "reassign": stage_reassign, "validate": stage_validate}
    if st == "figures":
        import ism_figures
        ism_figures.main()
    elif st == "all":
        for f in stages.values():
            f()
        import ism_figures
        ism_figures.main()
    else:
        stages[st]()
