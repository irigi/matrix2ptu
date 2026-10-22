"""Compare two converted acquisitions of the same sample area (e.g. the RAW and
the OBF acquisition of one dataset).  They are *separate* acquisitions, so this
measures agreement (rates, element fractions, decay shapes, image content after
registration), not identity.

Works on the detector-resolved PTU files and the conversion reports, i.e. on
the delivered outputs.

Usage:
    .venv/bin/python tools/compare_pair.py outputs/MatrixRAW/MatrixRAW outputs/MatrixOBF/MatrixOBF \
        [--json outputs/compare_Matrix.json]
(each argument is the output prefix: <prefix>_detector_resolved.ptu and
 <prefix>_conversion_report.json must exist)
"""
import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from matrix2ptu.validate import iter_ptu_records   # noqa: E402

N_EL = 23


def load(prefix):
    rep = json.load(open(prefix + "_conversion_report.json"))
    P = {k: v["value"] for k, v in rep["parameters"].items()}
    seqs = P["sequences"]
    unit = sum(s["line_accumulation"] for s in seqs)
    pattern = np.repeat(np.arange(len(seqs)), [s["line_accumulation"] for s in seqs])
    nx, ny = P["pixels_x"], P["pixels_y"]
    nseq = len(seqs)
    res_ps = P["microtime_resolution_ps"]
    nb = rep["derived"]["microtime_bins_per_period"]
    path = prefix + "_detector_resolved.ptu"

    # markers first (cheap second pass is simpler than buffering)
    ls, le = [], []
    for _, d, t in iter_ptu_records(path):
        m = d["marker"]
        ls.append(t[(m & 1) > 0])
        le.append(t[(m & 2) > 0])
    ls, le = np.concatenate(ls), np.concatenate(le)
    n_pass = ls.size

    img = np.zeros((nseq, ny, nx), np.int64)
    el_cnt = np.zeros((nseq, N_EL), np.int64)
    dec = np.zeros((nseq, nb), np.int64)
    dec_el = np.zeros((N_EL, nb), np.int64)          # sequence 0 only (confocal)
    pass_cnt = np.zeros(n_pass, np.int64)
    n_all = n_out = 0
    for _, d, t in iter_ptu_records(path):
        ph = (d["channel"] >= 0) & (d["marker"] == 0)
        tp, ch, dt = t[ph], d["channel"][ph].astype(np.int64), d["dtime"][ph].astype(np.int64)
        n_all += tp.size
        k = np.searchsorted(ls, tp, side="right") - 1
        kk = np.clip(k, 0, None)
        ok = (k >= 0) & (tp < le[kk])
        n_out += int(((k >= 0) & ~ok).sum())
        x = np.floor((tp - ls[kk]) / (le[kk] - ls[kk]) * nx).astype(np.int64)
        y = kk // unit
        ok &= (x >= 0) & (x < nx) & (y < ny)
        s = pattern[kk % unit]
        np.add.at(img, (s[ok], y[ok], x[ok]), 1)
        np.add.at(el_cnt, (s[ok], ch[ok]), 1)
        np.add.at(dec, (s[ok], np.minimum(dt[ok], nb - 1)), 1)
        m0 = ok & (s == 0)
        np.add.at(dec_el, (ch[m0], np.minimum(dt[m0], nb - 1)), 1)
        pass_cnt += np.bincount(kk[ok], minlength=n_pass)
    return dict(rep=rep, P=P, seqs=seqs, unit=unit, img=img, el_cnt=el_cnt, dec=dec,
                dec_el=dec_el, pass_cnt=pass_cnt, n_all=n_all, n_out=n_out, res_ps=res_ps,
                prefix=prefix, Tp_ns=rep["derived"]["laser_period_ns"])


def register(a, b):
    """Integer shift (dy, dx) that best aligns b to a (phase correlation), and
    the correlation of the overlapping parts after the shift."""
    A = np.fft.fft2(a - a.mean())
    B = np.fft.fft2(b - b.mean())
    R = A * np.conj(B)
    c = np.fft.ifft2(R / np.maximum(np.abs(R), 1e-12)).real
    dy, dx = np.unravel_index(c.argmax(), c.shape)
    dy = dy - a.shape[0] if dy > a.shape[0] // 2 else dy
    dx = dx - a.shape[1] if dx > a.shape[1] // 2 else dx
    ys = slice(max(dy, 0), a.shape[0] + min(dy, 0))
    xs = slice(max(dx, 0), a.shape[1] + min(dx, 0))
    yb = slice(max(-dy, 0), b.shape[0] + min(-dy, 0))
    xb = slice(max(-dx, 0), b.shape[1] + min(-dx, 0))
    r = float(np.corrcoef(a[ys, xs].ravel(), b[yb, xb].ravel())[0, 1])
    return int(dy), int(dx), r


def poisson_corr_ceiling(a, b):
    """Highest correlation two independent Poisson acquisitions of the same
    object can reach: sqrt(S_a/V_a * S_b/V_b), signal variance S = V - mean."""
    va, vb = a.var(), b.var()
    sa, sb = va - a.mean(), vb - b.mean()
    return float(np.sqrt(max(sa, 0) / va * max(sb, 0) / vb))


def decay_stats(h, res_ps, Tp_ns):
    h = h.astype(float)
    t = (np.arange(h.size) + 0.5) * res_ps * 1e-3     # ns
    w = 2 * np.pi / Tp_ns
    g = float((h * np.cos(w * t)).sum() / h.sum())
    s = float((h * np.sin(w * t)).sum() / h.sum())
    # smoothed peak and tail mean (from 0.5 ns after the peak, over 10 ns)
    k = 20
    hs = np.convolve(h, np.ones(k) / k, mode="same")
    p = int(hs.argmax())
    a0 = p + int(500 / res_ps)
    a1 = min(h.size, a0 + int(10000 / res_ps))
    tail_mean = float((h[a0:a1] * (t[a0:a1] - t[a0])).sum() / h[a0:a1].sum())
    return {"peak_ns": float(t[p]), "phasor_g": g, "phasor_s": s,
            "tail_mean_arrival_ns_0p5_to_10p5_after_peak": tail_mean}


def rising_edge_ps(h, res_ps):
    """Time (ps) where the smoothed decay first reaches half its maximum."""
    k = 10
    hs = np.convolve(h.astype(float), np.ones(k) / k, mode="same")
    p = int(hs.argmax())
    half = hs[p] / 2
    i = np.flatnonzero(hs[:p] < half)
    if not i.size:
        return float("nan")
    i = i[-1]
    frac = (half - hs[i]) / max(hs[i + 1] - hs[i], 1e-12)
    return float((i + frac + 0.5) * res_ps)


def to100(h, res_ps):
    f = int(round(100 / res_ps))
    n = h.size // f * f
    return h[:n].reshape(-1, f).sum(1)


def compare(a, b):
    out = {"a": a["prefix"], "b": b["prefix"]}
    Pa, Pb = a["P"], b["P"]
    out["timing"] = {k: [Pa.get(k), Pb.get(k)] for k in
                     ("laser_frequency_hz", "coarse_tick_ns", "sync_jitter_ps_rms",
                      "flyback_ms", "flyback_rate_ratio")}
    out["timing"]["missing_sync_records"] = [a["rep"]["derived"]["missing_sync_records_bridged"],
                                             b["rep"]["derived"]["missing_sync_records_bridged"]]
    out["photons_total"] = [a["n_all"], b["n_all"]]
    out["photons_outside_line_windows"] = [a["n_out"], b["n_out"]]
    seqs = []
    for i, s in enumerate(a["seqs"]):
        ia, ib = a["img"][i], b["img"][i]
        dy, dx, r_reg = register(ia.astype(float), ib.astype(float))
        r0 = float(np.corrcoef(ia.ravel(), ib.ravel())[0, 1])
        fa = a["el_cnt"][i] / a["el_cnt"][i].sum()
        fb = b["el_cnt"][i] / b["el_cnt"][i].sum()
        big = fa > 0.01
        da, db = a["dec"][i], b["dec"][i]
        A, B = to100(da, a["res_ps"]), to100(db, b["res_ps"])
        n = min(A.size, B.size)
        l1 = float(np.abs(A[:n] / A[:n].sum() - B[:n] / B[:n].sum()).sum())
        seqs.append({
            "name": s["name"], "line_accumulation": s["line_accumulation"],
            "image_photons": [int(ia.sum()), int(ib.sum())],
            "ratio_b_over_a": float(ib.sum() / ia.sum()),
            "image_corr_unregistered": r0,
            "registration_shift_px_dy_dx": [dy, dx],
            "image_corr_registered": r_reg,
            "image_corr_poisson_ceiling": poisson_corr_ceiling(ia.astype(float), ib.astype(float)),
            "element_fraction_max_abs_diff_pct_points": float(np.abs(fa - fb).max() * 100),
            "element_fraction_max_rel_diff_elements_gt1pct": float((np.abs(fa - fb) / fa)[big].max()),
            "element_fraction_corr": float(np.corrcoef(fa, fb)[0, 1]),
            "decay_L1_100ps": l1,
            "decay_stats": [decay_stats(da, a["res_ps"], a["Tp_ns"]),
                            decay_stats(db, b["res_ps"], b["Tp_ns"])],
        })
    out["sequences"] = seqs
    # per line pass: mean photons in each pass position of the acquisition cycle
    u = a["unit"]
    out["mean_photons_per_pass_position"] = [
        [float(x[k::u].mean()) for k in range(u)] for x in (a["pass_cnt"], b["pass_cnt"])]
    # per-element rising edge (sequence 0 = confocal), after delay correction
    ea = np.array([rising_edge_ps(h, a["res_ps"]) for h in a["dec_el"]])
    eb = np.array([rising_edge_ps(h, b["res_ps"]) for h in b["dec_el"]])
    w = a["dec_el"].sum(1)
    good = np.isfinite(ea) & np.isfinite(eb) & (w > 0)
    d = (eb - ea)[good]
    d0 = d - np.average(d, weights=w[good])
    out["element_rising_edge_ps"] = {
        "a": ea.round(1).tolist(), "b": eb.round(1).tolist(),
        "common_shift_b_minus_a_ps": float(np.average(d, weights=w[good])),
        "per_element_rms_after_common_shift_ps_weighted": float(np.sqrt(np.average(d0 ** 2, weights=w[good]))),
        "elements_used": int(good.sum())}
    # TDC tables (if both were self-calibrated)
    ta = a["prefix"] + "_tdc_calibration.npz"
    tb = b["prefix"] + "_tdc_calibration.npz"
    if os.path.exists(ta) and os.path.exists(tb):
        za, zb = np.load(ta), np.load(tb)
        ma, mb = str(za["method"]), str(zb["method"])
        tick_ps = Pa["coarse_tick_ns"] * 1e3
        rms = {}
        for k in za.files:
            if k.startswith("el") and k in zb.files:
                x, y = za[k], zb[k]
                both = ~np.isnan(x) & ~np.isnan(y)
                rms[k] = float(np.sqrt(np.mean(((x[both] - y[both]) * tick_ps) ** 2)))
        out["tdc_tables"] = {"methods": [ma, mb], "rms_diff_ps_per_channel": rms}
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("a")
    ap.add_argument("b")
    ap.add_argument("--json")
    args = ap.parse_args()
    A, B = load(args.a), load(args.b)
    out = compare(A, B)
    import re
    txt = json.dumps(out, indent=1)
    # keep number lists on one line
    txt = re.sub(r"\[\s*([^\[\]{}]*?)\s*\]", lambda m: "[" + re.sub(r"\s+", " ", m.group(1)) + "]", txt)
    if args.json:
        with open(args.json, "w") as f:
            f.write(txt)
    print(txt)


if __name__ == "__main__":
    main()
