"""Bead-based calibration of MATRIX detector-element shifts for pixel reassignment.

Input: per-element images (23, Y, X) of one acquisition sequence and,
for bead detection, the STED pooled image of the same field (STED separates
beads that merge in confocal).

Method (per sequence):
 1. detect beads as local maxima of the smoothed STED pooled image;
 2. classify each bead: isolated = no other detected bead (any brightness above
    the noise floor) within `isolation_nm`, window fully inside the image,
    single in STED (second-moment ellipticity below a limit);
 3. for each isolated bead and each element fit a symmetric 2-D Gaussian +
    constant background (Poisson maximum likelihood) in a square window;
    the bead reference position is the fit of the pooled image;
    element shift on that bead = element centre - pooled centre;
 4. the element shift = inverse-variance weighted mean over beads; its
    uncertainty is the larger of the statistical error and the
    bead-to-bead scatter / sqrt(n) (scatter-based errors are reported both ways).

Coordinates: x = column (scan direction), y = line (down), pixels of the
acquisition grid; shift s_e: element image ~ reference displaced by +s_e.
"""
from __future__ import annotations

import numpy as np
from scipy import ndimage, optimize

N_EL = 23
ROWS = [5, 4, 5, 4, 5]


def layout():
    """Imspector element layout (labels 1..23 = data element 0..22), units of pitch,
    x right, y down as drawn by Imspector."""
    pts = []
    for r, n in enumerate(ROWS):
        y = (r - 2) * np.sqrt(3) / 2
        for k in range(n):
            pts.append((k - (n - 1) / 2, y))
    return np.array(pts)


# --------------------------------------------------------------------------- detection
def detect_beads(img, sigma_px=1.5, min_sep_px=4, rel_thresh=0.12, noise_k=6.0):
    """Local maxima of the smoothed image.  Returns array (n, 3): y, x, peak value.
    A peak is kept if it exceeds rel_thresh x the 99.5th percentile and
    background + noise_k x background noise (robust estimates)."""
    sm = ndimage.gaussian_filter(img.astype(float), sigma_px)
    mx = ndimage.maximum_filter(sm, size=2 * min_sep_px + 1)
    bg = np.median(sm)
    noise = 1.4826 * np.median(np.abs(sm - bg))
    thr = max(rel_thresh * np.percentile(sm, 99.5), bg + noise_k * max(noise, 1e-9))
    yy, xx = np.nonzero((sm == mx) & (sm > thr))
    return np.c_[yy, xx, sm[yy, xx]], {"threshold": float(thr), "background": float(bg), "noise": float(noise)}


def second_moments(win):
    """Background-subtracted second-moment ellipticity (1 - minor/major)."""
    w = win - np.percentile(win, 20)
    w[w < 0] = 0
    if w.sum() <= 0:
        return np.nan
    yy, xx = np.indices(w.shape)
    m = w.sum()
    cy, cx = (w * yy).sum() / m, (w * xx).sum() / m
    cyy = (w * (yy - cy) ** 2).sum() / m
    cxx = (w * (xx - cx) ** 2).sum() / m
    cxy = (w * (yy - cy) * (xx - cx)).sum() / m
    ev = np.linalg.eigvalsh([[cxx, cxy], [cxy, cyy]])
    return float(1 - np.sqrt(max(ev[0], 1e-12) / ev[1]))


def classify_beads(peaks, shape, isolation_px, half_win_px, sted_img=None, max_ellipticity=0.35,
                   sted_half_win=5, all_sources=None, min_peak=0.0):
    """Return a list of dicts (one per peak) with position, nearest-neighbour
    distance (to any source in `all_sources`, default: the peaks themselves) and
    flags.  isolated = nn >= isolation_px, window inside the image, single in
    STED (ellipticity <= max_ellipticity) and peak >= min_peak."""
    H, W = shape
    src = (peaks if all_sources is None else all_sources)[:, :2].astype(float)
    out = []
    for i, (y, x, v) in enumerate(peaks):
        d = np.hypot(src[:, 0] - y, src[:, 1] - x)
        d = d[d > 0.5]
        nn = float(d.min()) if d.size else np.inf
        edge = (y - half_win_px < 0 or x - half_win_px < 0 or y + half_win_px >= H or x + half_win_px >= W)
        ell = np.nan
        if sted_img is not None:
            y0, x0 = int(y), int(x)
            win = sted_img[max(0, y0 - sted_half_win):y0 + sted_half_win + 1,
                           max(0, x0 - sted_half_win):x0 + sted_half_win + 1].astype(float)
            ell = second_moments(win)
        single = not (ell > max_ellipticity)
        out.append({"y": float(y), "x": float(x), "peak": float(v), "nn_px": nn, "edge": bool(edge),
                    "sted_ellipticity": ell, "single_in_sted": bool(single),
                    "bright": bool(v >= min_peak),
                    "isolated": bool(nn >= isolation_px and not edge and single and v >= min_peak)})
    return out


# --------------------------------------------------------------------------- fitting
def _gauss(p, yy, xx):
    a, y0, x0, s, b = p
    return a * np.exp(-((yy - y0) ** 2 + (xx - x0) ** 2) / (2 * s * s)) + b


def fit_gauss(win, y_guess, x_guess, s_guess, neighbours=()):
    """Poisson ML fit of symmetric Gaussians (common sigma) + constant background.
    The first Gaussian is the bead of interest; `neighbours` are (y, x) start
    positions (window coordinates) of other beads whose light reaches the window;
    they are fitted jointly (own amplitude and position) so they do not bias the
    bead.  Returns (params [amp, y, x, sigma, bg] of the bead, their 1-sigma
    errors from the Poisson Fisher information, success flag)."""
    yy, xx = np.indices(win.shape, dtype=float)
    d = win.astype(float).ravel()
    yy, xx = yy.ravel(), xx.ravel()
    H, W = win.shape
    src = [(y_guess, x_guess)] + [tuple(n) for n in neighbours]
    ns = len(src)
    b0 = max(np.percentile(d, 10), 0.01)
    a0 = max(d.max() - b0, 1.0)
    # p = [sigma, bg, a_1, y_1, x_1, a_2, ...]
    p0 = [s_guess, b0]
    for (y, x) in src:
        yi = int(np.clip(round(y), 0, H - 1)); xi = int(np.clip(round(x), 0, W - 1))
        p0 += [max(win[yi, xi] - b0, 0.1 * a0), y, x]
    p0 = np.array(p0, float)

    def parts(p):
        s, b = p[0], p[1]
        A = p[2::3]; Y = p[3::3]; X = p[4::3]
        r2 = (yy[None] - Y[:, None]) ** 2 + (xx[None] - X[:, None]) ** 2
        g = np.exp(-r2 / (2 * s * s))
        return s, b, A, Y, X, r2, g

    def model(p):
        s, b, A, Y, X, r2, g = parts(p)
        return (A[:, None] * g).sum(0) + b

    def grad_model(p):
        s, b, A, Y, X, r2, g = parts(p)
        rows = [(A[:, None] * g * r2 / s ** 3).sum(0), np.ones_like(yy)]
        for i in range(ns):
            rows += [g[i], A[i] * g[i] * (yy - Y[i]) / s ** 2, A[i] * g[i] * (xx - X[i]) / s ** 2]
        return np.stack(rows)

    def nll(p):
        m = model(p)
        if (m <= 0).any():
            return 1e30
        return float((m - d * np.log(m)).sum())

    def jac(p):
        m = model(p)
        return (grad_model(p) * (1 - d / np.maximum(m, 1e-12))).sum(1)

    bounds = [(0.4, max(H, W)), (1e-4, None)]
    for (y, x) in src:
        bounds += [(1e-3, None), (-3, H + 2), (-3, W + 2)]
    r = optimize.minimize(nll, p0, jac=jac, method="L-BFGS-B", bounds=bounds,
                          options={"maxiter": 2000})
    p = r.x
    m = model(p)
    G = grad_model(p)
    F = (G / np.maximum(m, 1e-12)) @ G.T          # Fisher information (Poisson)
    try:
        cov = np.linalg.inv(F)
        e = np.sqrt(np.abs(np.diag(cov)))
    except np.linalg.LinAlgError:
        e = np.full(p.size, np.nan)
    out = np.array([p[2], p[3], p[4], p[0], p[1]])
    err = np.array([e[2], e[3], e[4], e[0], e[1]])
    ok = bool(r.success) and 0 <= p[3] <= H - 1 and 0 <= p[4] <= W - 1
    return out, err, ok


def centroid(win, bg):
    w = win.astype(float) - bg
    yy, xx = np.indices(w.shape)
    m = w.sum()
    return (w * yy).sum() / m, (w * xx).sum() / m


def bead_element_shifts(stack, beads, half_win_px, s_guess_px, all_sources=None, reach_px=None):
    """Per bead and element: Gaussian-fit centre minus pooled-image centre.
    Sources from `all_sources` (y, x) within the window + `reach_px` are fitted
    jointly as neighbours.  Returns one dict per bead with arrays (23, 2) [x, y]:
    shift, err, centroid_shift (centroid of the background-subtracted window;
    unreliable when neighbours are present), element signal and sigma."""
    pooled = stack.sum(0)
    reach = half_win_px + (reach_px if reach_px is not None else 2 * s_guess_px)
    rows = []
    for b in beads:
        y0, x0 = int(round(b["y"])), int(round(b["x"]))
        oy, ox = y0 - half_win_px, x0 - half_win_px
        sl = (slice(oy, y0 + half_win_px + 1), slice(ox, x0 + half_win_px + 1))
        nb = []
        if all_sources is not None:
            for (yy, xx) in np.asarray(all_sources)[:, :2]:
                if 0.5 < np.hypot(yy - b["y"], xx - b["x"]) and max(abs(yy - y0), abs(xx - x0)) <= reach:
                    nb.append((yy - oy, xx - ox))
        pw = pooled[sl]
        pp, pe, ok = fit_gauss(pw, b["y"] - oy, b["x"] - ox, s_guess_px, nb)
        if not ok or not np.isfinite(pe).all():
            continue
        cy_p, cx_p = centroid(pw, pp[4])
        sh = np.zeros((N_EL, 2)); er = np.zeros((N_EL, 2)); cs = np.zeros((N_EL, 2))
        n_el = np.zeros(N_EL); sig = np.zeros(N_EL); good = np.ones(N_EL, bool)
        for e in range(N_EL):
            w = stack[e][sl]
            p, err, ok_e = fit_gauss(w, pp[1], pp[2], pp[3], nb)
            good[e] = ok_e and np.isfinite(err).all()
            sh[e] = (p[2] - pp[2], p[1] - pp[1])
            er[e] = (err[2], err[1])
            cy, cx = centroid(w, p[4])
            cs[e] = (cx - cx_p, cy - cy_p)
            n_el[e] = 2 * np.pi * p[0] * p[3] ** 2
            sig[e] = p[3]
        rows.append({"bead": b, "n_neighbours_in_fit": len(nb), "pooled_fit": pp, "pooled_err": pe,
                     "shift": sh, "err": er, "fit_ok": good,
                     "centroid_shift": cs, "element_signal": n_el, "element_sigma": sig,
                     "window_origin": (oy, ox)})
    return rows


def combine(rows, use="shift"):
    """Inverse-variance weighted mean over beads (element fits that failed get
    zero weight).  Returns mean (23, 2), statistical error, bead-to-bead
    scatter error (std / sqrt(n)), reduced chi2 of the scatter."""
    S = np.array([r[use] for r in rows])                # (n, 23, 2)
    E = np.array([r["err"] for r in rows])
    ok = np.array([r["fit_ok"] for r in rows])[..., None] & np.isfinite(E) & np.isfinite(S)
    w = np.where(ok, 1 / np.maximum(np.where(ok, E, 1), 1e-6) ** 2, 0.0)
    S = np.where(ok, S, 0.0)
    mean = (w * S).sum(0) / w.sum(0)
    stat = 1 / np.sqrt(w.sum(0))
    n = ok.sum(0)
    dev = np.where(ok, S - mean, 0.0)
    scat = np.sqrt((dev ** 2).sum(0) / np.maximum(n - 1, 1)) / np.sqrt(np.maximum(n, 1))
    chi2 = (dev ** 2 * w).sum(0) / np.maximum(n - 1, 1)
    return mean, stat, scat, chi2


def per_bead_scale(rows, ref):
    """Scale of each bead's shift pattern relative to `ref` (affine fit)."""
    return np.array([describe_affine(fit_affine(ref, r["shift"], np.ones(N_EL))[0])["scale_per_pitch"]
                     for r in rows])


# --------------------------------------------------------------------------- lattice model
def fit_affine(P, S, w):
    X = np.c_[P, np.ones(len(P))]
    W = np.sqrt(w)[:, None]
    coef, *_ = np.linalg.lstsq(X * W, S * W, rcond=None)
    A, t = coef[:2].T, coef[2]
    return A, t, S - (P @ A.T + t)


def describe_affine(A):
    U, sv, Vt = np.linalg.svd(A)
    mirror = bool(np.linalg.det(A) < 0)
    R = U @ (np.diag([1, -1]) if mirror else np.eye(2)) @ Vt
    return {"scale_per_pitch": float(np.sqrt(abs(np.linalg.det(A)))),
            "singular_values": sv.tolist(), "aspect_min_over_max": float(sv[1] / sv[0]),
            "rotation_deg_x_right_y_down": float(np.degrees(np.arctan2(R[1, 0], R[0, 0])) % 360),
            "mirrored": mirror}


# --------------------------------------------------------------------------- image-domain tools
def fourier_shift(img, sx, sy):
    """Shift image content by (sx, sy) pixels (sub-pixel, periodic)."""
    F = np.fft.fft2(img)
    ky = np.fft.fftfreq(img.shape[0])[:, None]
    kx = np.fft.fftfreq(img.shape[1])[None, :]
    return np.real(np.fft.ifft2(F * np.exp(-2j * np.pi * (ky * sy + kx * sx))))


def xcorr_shift(img, ref, win=10):
    """Shift (sx, sy) with img(r) ~ ref(r - s), from the peak of the circular
    cross-correlation (3-point Gaussian sub-pixel interpolation).

    Circular, not zero-padded: with zero padding, the mean-subtracted images
    produce a triangular pedestal peaked at zero shift that biases the estimate
    towards zero (strongly for sparse images).  Wrap-around affects only a
    border as wide as the shift (a few pixels)."""
    H, W = img.shape
    a = img - img.mean()
    b = ref - ref.mean()
    c = np.fft.fftshift(np.fft.irfft2(np.fft.rfft2(a) * np.conj(np.fft.rfft2(b)), s=(H, W)))
    cy, cx = H // 2, W // 2
    sub = c[cy - win:cy + win + 1, cx - win:cx + win + 1]
    i, j = np.unravel_index(sub.argmax(), sub.shape)
    i, j = i + cy - win, j + cx - win

    def sp(m, o, p):
        if min(m, o, p) > 0:
            m, o, p = np.log(m), np.log(o), np.log(p)
        den = m - 2 * o + p
        return 0.0 if den == 0 else 0.5 * (m - p) / den
    return j + sp(c[i, j - 1], c[i, j], c[i, j + 1]) - cx, i + sp(c[i - 1, j], c[i, j], c[i + 1, j]) - cy


def frc(a, b, n_rings=None):
    """Fourier ring correlation of two independent images; returns (freq [1/px], frc)."""
    H, W = a.shape
    wy = np.hanning(H)[:, None] * np.hanning(W)[None, :]
    A = np.fft.fftshift(np.fft.fft2((a - a.mean()) * wy))
    B = np.fft.fftshift(np.fft.fft2((b - b.mean()) * wy))
    yy, xx = np.indices((H, W))
    r = np.hypot(yy - H // 2, xx - W // 2).astype(int)
    n = n_rings or min(H, W) // 2
    num = np.bincount(r.ravel(), (A * np.conj(B)).real.ravel(), minlength=n)[:n]
    da = np.bincount(r.ravel(), (np.abs(A) ** 2).ravel(), minlength=n)[:n]
    db = np.bincount(r.ravel(), (np.abs(B) ** 2).ravel(), minlength=n)[:n]
    return np.arange(n) / min(H, W), num / np.sqrt(da * db + 1e-30)


def frc_resolution(freq, curve, thr=1 / 7, smooth=3):
    c = np.convolve(curve, np.ones(smooth) / smooth, mode="same")
    for i in range(3, len(c)):
        if c[i] < thr:
            f = freq[i - 1] + (freq[i] - freq[i - 1]) * (c[i - 1] - thr) / (c[i - 1] - c[i])
            return 1 / f
    return np.nan


# --------------------------------------------------------------------------- calibration driver
def _crop(img, y0, x0, h):
    return img[..., y0 - h:y0 + h + 1, x0 - h:x0 + h + 1]


def loo_validation(stack, rows, all_sources, half_win_px, s_guess_px, scales=(1.0,)):
    """Leave-one-bead-out test: for each bead, shifts from all OTHER beads (times
    each factor in `scales`) are applied to that bead's element images (Fourier
    sub-pixel shift) and the pooled bead is refitted.  Returns per-bead FWHM and
    peak ratios after/before, and the central-element FWHM for reference."""
    hw, pad = half_win_px, half_win_px + 8
    out = []
    for i, r in enumerate(rows):
        m = combine([q for j, q in enumerate(rows) if j != i])[0]
        b = r["bead"]
        y0, x0 = int(round(b["y"])), int(round(b["x"]))
        sub = _crop(stack, y0, x0, pad)
        if sub.shape[1:] != (2 * pad + 1, 2 * pad + 1):
            continue
        nb = [(yy - y0 + hw, xx - x0 + hw) for yy, xx, _ in all_sources
              if 0.5 < np.hypot(yy - b["y"], xx - b["x"]) and max(abs(yy - y0), abs(xx - x0)) <= hw + 2 * s_guess_px]
        c = slice(pad - hw, pad + hw + 1)
        yg, xg = b["y"] - y0 + hw, b["x"] - x0 + hw
        p0 = fit_gauss(sub.sum(0)[c, c], yg, xg, s_guess_px, nb)[0]
        pc = fit_gauss(sub[11][c, c], yg, xg, s_guess_px, nb)[0]
        rec = {"x": b["x"], "y": b["y"], "sigma_before_px": float(p0[3]), "sigma_central_element_px": float(pc[3])}
        for k in scales:
            a = sum(fourier_shift(sub[e], -k * m[e, 0], -k * m[e, 1]) for e in range(N_EL))
            pa = fit_gauss(a[c, c], yg, xg, s_guess_px, nb)[0]
            rec[f"fwhm_ratio_k{k:g}"] = float(pa[3] / p0[3])
            rec[f"peak_ratio_k{k:g}"] = float(pa[0] / p0[0])
        out.append(rec)
    return out


def permutation_test(P, S, w, n=2000, seed=0):
    rng = np.random.default_rng(seed)
    ref = float(np.sqrt(np.average((fit_affine(P, S, w)[2] ** 2).sum(1), weights=w)))
    vals = []
    for _ in range(n):
        pm = rng.permutation(N_EL)
        vals.append(np.sqrt(np.average((fit_affine(P[pm], S, w)[2] ** 2).sum(1), weights=w)))
    vals = np.array(vals)
    return {"wrms_identity_mapping_px": ref, "random_permutations": n,
            "fraction_of_permutations_fitting_as_well": float((vals <= ref).mean()),
            "median_wrms_random_px": float(np.median(vals))}


def calibrate_sequence(stack, sted_pooled, pixel_size_nm, half_win_px, s_guess_px,
                       isolation_px=16, max_neighbours=1, robustness_radii=(12, 20),
                       mad_k=4.0, loo_scales=(0.6, 0.8, 1.0, 1.2, 1.4, 1.6)):
    """Full bead calibration of one acquisition sequence.  See module docstring."""
    all_src, det_info = detect_beads(sted_pooled, rel_thresh=0.0, noise_k=5)
    cand, _ = detect_beads(sted_pooled, rel_thresh=0.05)
    edge = half_win_px + 8

    def select(R):
        cl = classify_beads(cand, sted_pooled.shape, R, edge, sted_pooled, all_sources=all_src)
        iso = [c for c in cl if c["isolated"]]
        rows = bead_element_shifts(stack, iso, half_win_px, s_guess_px, all_sources=all_src)
        crowded = [r for r in rows if r["n_neighbours_in_fit"] > max_neighbours]
        rows = [r for r in rows if r["n_neighbours_in_fit"] <= max_neighbours]
        return cl, rows, crowded

    table, rows, crowded = select(isolation_px)
    # robust outlier rejection on the per-bead scale of the shift pattern
    prov = combine(rows)[0]
    k = per_bead_scale(rows, prov)
    med = np.median(k)
    mad = 1.4826 * np.median(np.abs(k - med))
    bad = np.abs(k - med) > mad_k * max(mad, 0.02)
    rejected = [{"x": r["bead"]["x"], "y": r["bead"]["y"], "scale": float(kk)} for r, kk, b in zip(rows, k, bad) if b]
    rows = [r for r, b in zip(rows, bad) if not b]
    mean, stat, scat, chi2 = combine(rows)
    unc = np.maximum(stat, scat)

    # robustness: isolation radius
    rob = {}
    for R in robustness_radii:
        _, rr, _ = select(R)
        if len(rr) >= 3:
            m2 = combine(rr)[0]
            rob[f"R{R}px"] = {"n_beads": len(rr), "rms_diff_px": float(np.sqrt(((m2 - mean) ** 2).sum(1).mean()))}
    # split halves: left/right of the field and alternating beads
    halves = {}
    xs = np.array([r["bead"]["x"] for r in rows])
    for name, sel in (("left_vs_right", xs < np.median(xs)), ("alternate", np.arange(len(rows)) % 2 == 0)):
        A = [r for r, s_ in zip(rows, sel) if s_]
        B = [r for r, s_ in zip(rows, sel) if not s_]
        if len(A) >= 2 and len(B) >= 2:
            ma, ea = combine(A)[0], np.maximum(*combine(A)[1:3])
            mb, eb = combine(B)[0], np.maximum(*combine(B)[1:3])
            z = (ma - mb) / np.sqrt(ea ** 2 + eb ** 2)
            halves[name] = {"n": [len(A), len(B)], "rms_diff_px": float(np.sqrt(((ma - mb) ** 2).sum(1).mean())),
                            "rms_z": float(np.sqrt((z ** 2).mean())),
                            "scale_B_over_A": describe_affine(fit_affine(ma, mb, np.ones(N_EL))[0])["scale_per_pitch"]}
    # field dependence: per-bead scale vs position
    kk = per_bead_scale(rows, mean)
    X = np.c_[np.ones(len(rows)), [r["bead"]["x"] for r in rows], [r["bead"]["y"] for r in rows]]
    coef, *_ = np.linalg.lstsq(X, kk, rcond=None)
    # lattice diagnostics
    P = layout()
    wts = np.array([r["element_signal"] for r in rows]).sum(0)
    A, t, res = fit_affine(P, mean, wts)
    loo = loo_validation(stack, rows, all_src, half_win_px, s_guess_px, loo_scales)

    def summ(key):
        v = np.array([l[key] for l in loo])
        return {"median": float(np.median(v)), "mean": float(v.mean()),
                "sem": float(v.std(ddof=1) / np.sqrt(len(v))) if len(v) > 1 else None,
                "n_improved": int((v < 1).sum()) if "fwhm" in key else int((v > 1).sum()), "n": len(v)}
    return {
        "pixel_size_nm": pixel_size_nm,
        "method": {"detection": "local maxima of the smoothed STED pooled image", "isolation_px": isolation_px,
                   "max_neighbours_in_fit": max_neighbours, "fit_half_window_px": half_win_px,
                   "fit": "Poisson ML symmetric 2-D Gaussian + background, neighbours fitted jointly",
                   "combination": "inverse-variance weighted mean over beads",
                   "outlier_rejection": f"per-bead pattern scale beyond {mad_k} x MAD"},
        "detection": det_info, "n_sources_detected": int(len(all_src)), "n_candidates": int(len(cand)),
        "bead_table": table,
        "beads_used": [{"x": r["bead"]["x"], "y": r["bead"]["y"], "n_neighbours_in_fit": r["n_neighbours_in_fit"],
                        "pooled_sigma_px": float(r["pooled_fit"][3]), "pooled_amplitude": float(r["pooled_fit"][0]),
                        "pattern_scale": float(s_)} for r, s_ in zip(rows, kk)],
        "beads_excluded_crowded": [{"x": r["bead"]["x"], "y": r["bead"]["y"], "n_neighbours_in_fit": r["n_neighbours_in_fit"]} for r in crowded],
        "beads_rejected_outlier": rejected,
        "shift_px": mean, "uncertainty_px": unc, "stat_err_px": stat, "scatter_err_px": scat, "chi2_scatter": chi2,
        "per_bead_shift_px": np.array([r["shift"] for r in rows]),
        "element_photons_in_beads": wts,
        "robustness_isolation_radius": rob, "split_halves": halves,
        "field_dependence": {"model": "pattern scale = c0 + cx * x + cy * y (x, y in px)",
                             "coef": coef.tolist(), "scale_change_across_field": float(abs(coef[1]) * 250 + abs(coef[2]) * 250),
                             "per_bead_scale_sd": float(kk.std(ddof=1))},
        "lattice_fit": {**describe_affine(A), "translation_px": t.tolist(),
                        "residual_rms_px": float(np.sqrt((res ** 2).sum(1).mean())),
                        "uncertainty_rms_px": float(np.sqrt((unc ** 2).sum(1).mean())),
                        "label_mapping_test": permutation_test(P, mean, wts)},
        "loo_validation": {"per_bead": loo,
                           "summary": {f"k{k:g}": {"fwhm_ratio": summ(f"fwhm_ratio_k{k:g}"),
                                                   "peak_ratio": summ(f"peak_ratio_k{k:g}")} for k in loo_scales},
                           "fwhm_before_nm_median": float(np.median([l["sigma_before_px"] for l in loo]) * 2.3548 * pixel_size_nm),
                           "fwhm_central_element_nm_median": float(np.median([l["sigma_central_element_px"] for l in loo]) * 2.3548 * pixel_size_nm)},
    }
