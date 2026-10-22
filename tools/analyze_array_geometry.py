"""Assessment of MATRIX detector-element geometry for future pixel reassignment.

Analysis only - no reassignment is performed.  For each dataset, per-element
intensity images are formed and the image shift of every element relative to
the sum of all *other* elements is measured (cross-correlation, sub-pixel peak).
These "fingerprint" shifts are what adaptive pixel reassignment (APR) uses.
The shift pattern is then fitted with the hexagonal 23-element layout shown by
Imspector (rows of 5-4-5-4-5, labels 1..23 left->right, top->bottom; hypothesis:
data element = label - 1):

    shift_i = A @ p_i + t        (p_i: layout position in units of the pitch)

A general 2x2 A (affine: scale, rotation, aspect, shear, mirror) and a
similarity A (scale, rotation, optional mirror) are both fitted, weighted by
the element photon counts.  The Huygens parameters (spacing, aspect ratio,
rotation, mirroring) are compared with the fitted lattice under several
candidate conventions.

Usage:
    .venv/bin/python tools/analyze_array_geometry.py [--json outputs/array_geometry.json]
        [--fig docs/figures/f7_array_shifts.png]
"""
import argparse
import itertools
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from matrix2ptu.obf import read_obf                 # noqa: E402
from matrix2ptu.validate import iter_ptu_records    # noqa: E402

ROOT = os.path.join(os.path.dirname(__file__), "..")
N_EL = 23
ROWS = [5, 4, 5, 4, 5]
HUYGENS = {"spacing_nm": 108.1, "aspect_x_to_y": 0.7918, "rotation_deg": 274.505, "mirroring": "none"}


def layout():
    """Imspector layout, units of the element pitch; x right, y down (as drawn)."""
    pts = []
    for r, n in enumerate(ROWS):
        y = (r - 2) * np.sqrt(3) / 2
        for k in range(n):
            pts.append((k - (n - 1) / 2, y))
    return np.array(pts)


# ----------------------------------------------------------------------------- images
def images_from_lightbox(path):
    st, _, _ = read_obf(path)
    s = st[0]
    a = s.read_data()[0, :, 0].astype(float)        # (23, y, x)
    px_nm = s.lengths_m[0] / a.shape[2] * 1e9
    return {"all": a}, px_nm


def images_from_ptu(prefix):
    """Per-element images (sequence, element, y, x) from a detector-resolved PTU."""
    rep = json.load(open(prefix + "_conversion_report.json"))
    P = {k: v["value"] for k, v in rep["parameters"].items()}
    seqs = P["sequences"]
    unit = sum(s["line_accumulation"] for s in seqs)
    pattern = np.repeat(np.arange(len(seqs)), [s["line_accumulation"] for s in seqs])
    nx, ny = P["pixels_x"], P["pixels_y"]
    path = prefix + "_detector_resolved.ptu"
    ls, le = [], []
    for _, d, t in iter_ptu_records(path):
        ls.append(t[(d["marker"] & 1) > 0])
        le.append(t[(d["marker"] & 2) > 0])
    ls, le = np.concatenate(ls), np.concatenate(le)
    img = np.zeros((len(seqs), N_EL, ny, nx), np.int32)
    for _, d, t in iter_ptu_records(path):
        ph = (d["channel"] >= 0) & (d["marker"] == 0)
        tp, ch = t[ph], d["channel"][ph].astype(np.int64)
        k = np.searchsorted(ls, tp, side="right") - 1
        kk = np.clip(k, 0, None)
        ok = (k >= 0) & (tp < le[kk])
        x = np.floor((tp - ls[kk]) / (le[kk] - ls[kk]) * nx).astype(np.int64)
        y = kk // unit
        ok &= (x >= 0) & (x < nx) & (y < ny)
        np.add.at(img, (pattern[kk % unit][ok], ch[ok], y[ok], x[ok]), 1)
    return {s["name"]: img[i].astype(float) for i, s in enumerate(seqs)}, P["pixel_size_um"] * 1e3


# ----------------------------------------------------------------------------- shifts
def _peak_subpixel(c, i, j):
    def fit(m, o, p):
        if min(m, o, p) > 0:                         # Gaussian 3-point
            lm, lo, lp = np.log(m), np.log(o), np.log(p)
            den = lm - 2 * lo + lp
        else:
            lm, lo, lp, den = m, o, p, m - 2 * o + p
        return 0.0 if den == 0 else 0.5 * (lm - lp) / den
    return i + fit(c[i - 1, j], c[i, j], c[i + 1, j]), j + fit(c[i, j - 1], c[i, j], c[i, j + 1])


def xcorr_shift(img, ref):
    """Shift (dy, dx) in pixels such that img(r) ~ ref(r - shift).  Linear
    (zero-padded) cross-correlation of mean-subtracted images."""
    H, W = img.shape
    a = img - img.mean()
    b = ref - ref.mean()
    F = np.fft.rfft2(a, s=(2 * H, 2 * W)) * np.conj(np.fft.rfft2(b, s=(2 * H, 2 * W)))
    c = np.fft.fftshift(np.fft.irfft2(F, s=(2 * H, 2 * W)))
    c0 = c.copy()
    cy, cx = H, W
    win = 12                                          # search +-12 px around zero shift
    sub = c0[cy - win:cy + win + 1, cx - win:cx + win + 1]
    i, j = np.unravel_index(sub.argmax(), sub.shape)
    i, j = i + cy - win, j + cx - win
    if c[i, j] <= 0:
        c = c - c.min() + 1e-9
    y, x = _peak_subpixel(c, i, j)
    return y - cy, x - cx


def element_shifts(stack, px_nm):
    total = stack.sum(0)
    out = np.zeros((N_EL, 2))
    for e in range(N_EL):
        dy, dx = xcorr_shift(stack[e], total - stack[e])   # leave-one-out reference
        out[e] = (dx * px_nm, dy * px_nm)                   # (x, y) in nm, y = line index (down)
    return out, stack.reshape(N_EL, -1).sum(1)


# ----------------------------------------------------------------------------- fits
def fit_affine(P, S, w):
    X = np.c_[P, np.ones(len(P))]
    W = np.sqrt(w)[:, None]
    coef, *_ = np.linalg.lstsq(X * W, S * W, rcond=None)
    A, t = coef[:2].T, coef[2]
    res = S - (P @ A.T + t)
    return A, t, res


def fit_similarity(P, S, w, allow_mirror=True):
    """Weighted Procrustes: S ~ s R P (+ mirror) + t."""
    wn = w / w.sum()
    mp, ms = wn @ P, wn @ S
    Pc, Sc = P - mp, S - ms
    C = (Sc * wn[:, None]).T @ Pc
    U, sv, Vt = np.linalg.svd(C)
    d = np.sign(np.linalg.det(U @ Vt))
    if not allow_mirror:
        D = np.diag([1, d])
    else:
        D = np.eye(2)
    R = U @ D @ Vt
    scale = np.trace(np.diag(sv) @ D) / (wn @ (Pc ** 2).sum(1))
    A = scale * R
    t = ms - A @ mp
    res = S - (P @ A.T + t)
    return A, t, res


def describe(A):
    """Scale (geometric mean), rotation of the best-fitting rotation, mirror flag,
    singular-value ratio (aspect) and its axis angle."""
    U, sv, Vt = np.linalg.svd(A)
    mirror = bool(np.linalg.det(A) < 0)
    D = np.diag([1, -1]) if mirror else np.eye(2)
    R = U @ D @ Vt if mirror else U @ Vt
    rot = np.degrees(np.arctan2(R[1, 0], R[0, 0])) % 360
    return {"scale_nm_per_pitch": float(np.sqrt(abs(np.linalg.det(A)))),
            "singular_values_nm": sv.round(2).tolist(),
            "aspect_min_over_max": float(sv[1] / sv[0]),
            "rotation_deg_(x_right_y_down_clockwise_positive)": float(rot),
            "mirrored": mirror}


def wrms(res, w):
    return float(np.sqrt(np.average((res ** 2).sum(1), weights=w)))


def permutation_test(P, S, w, n=2000, seed=0):
    rng = np.random.default_rng(seed)
    vals = []
    for _ in range(n):
        p = rng.permutation(N_EL)
        _, _, r = fit_affine(P[p], S, w)
        vals.append(wrms(r, w))
    return np.array(vals)


def huygens_conventions(P, S, w):
    """Model positions from the Huygens numbers under candidate conventions;
    fit only a common scale k and offset (reassignment factor), report residual."""
    th = np.radians(HUYGENS["rotation_deg"])
    a = HUYGENS["aspect_x_to_y"]
    rows = []
    for sign, ax, flip, swap in itertools.product((1, -1), ("x", "y"), (False, True), (False, True)):
        Q = P[:, ::-1] if swap else P.copy()          # swap: layout rotated by 90 deg (rows <-> columns)
        Dm = np.diag([a, 1.0]) if ax == "x" else np.diag([1.0, a])
        R = np.array([[np.cos(sign * th), -np.sin(sign * th)], [np.sin(sign * th), np.cos(sign * th)]])
        M = R @ Dm
        Hp = Q @ M.T * HUYGENS["spacing_nm"]
        if flip:
            Hp[:, 1] *= -1
        wn = w / w.sum()
        Hc, Sc = Hp - wn @ Hp, S - wn @ S
        k = float((wn * (Hc * Sc).sum(1)).sum() / (wn * (Hc ** 2).sum(1)).sum())
        res = Sc - k * Hc
        rows.append({"rotation_sign": sign, "aspect_on": ax, "y_flip": flip, "layout_rows_vertical": swap,
                     "factor_k": k, "residual_nm_wrms": wrms(res, w)})
    rows.sort(key=lambda r: r["residual_nm_wrms"])
    return rows


def intensity_centroid(P, counts):
    f = counts / counts.sum()
    return (f @ P).tolist(), f.round(5).tolist()


# ----------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", default=os.path.join(ROOT, "outputs", "array_geometry.json"))
    ap.add_argument("--fig", default=os.path.join(ROOT, "docs", "figures", "f7_array_shifts.png"))
    args = ap.parse_args()
    P = layout()
    sets = []
    lb, px = images_from_lightbox(os.path.join(ROOT, "inputs2", "MatrixLB.obf"))
    sets.append(("MatrixLB.obf (Lightbox, 561 nm, confocal)", lb["all"], px))
    for prefix, label in (("OrangebeadsRAW/Sequence_orangebeadsRAW", "Orangebeads RAW"),
                          ("OrangebeadsOBF/Sequence_orangebeadsOBF", "Orangebeads OBF"),
                          ("MatrixRAW/MatrixRAW", "Matrix 640 RAW"),
                          ("MatrixOBF/MatrixOBF", "Matrix 640 OBF")):
        ims, px = images_from_ptu(os.path.join(ROOT, "outputs", prefix))
        for name, st in ims.items():
            sets.append((f"{label} {name}", st, px))

    out = {"layout": "Imspector hex_23, rows 5-4-5-4-5, element = label - 1 (tested below)",
           "huygens_reported": HUYGENS, "datasets": []}
    shifts_all = {}
    for label, stack, px in sets:
        S, counts = element_shifts(stack, px)
        w = counts.astype(float)
        A, t, res = fit_affine(P, S, w)
        As, ts, ress = fit_similarity(P, S, w)
        perm = permutation_test(P, S, w)
        d = {"dataset": label, "pixel_nm": px, "photons": int(counts.sum()),
             "photons_per_element": counts.astype(int).tolist(),
             "shifts_nm_xy": S.round(1).tolist(),
             "affine": describe(A) | {"residual_nm_wrms": wrms(res, w)},
             "similarity": describe(As) | {"residual_nm_wrms": wrms(ress, w)},
             "layout_mapping_test": {
                 "residual_label_minus_1_nm": wrms(res, w),
                 "residual_random_permutations_nm_min_median": [float(perm.min()), float(np.median(perm))],
                 "fraction_of_random_permutations_fitting_as_well": float((perm <= wrms(res, w)).mean())},
             "huygens_conventions_best3": huygens_conventions(P, S, w)[:3],
             "intensity_centroid_pitch_units_xy": intensity_centroid(P, w)[0]}
        out["datasets"].append(d)
        shifts_all[label] = (S, w, A, t)
        print(f"{label:52s} px {px:4.0f} nm  photons {int(counts.sum()):9d}  "
              f"affine scale {d['affine']['scale_nm_per_pitch']:5.1f} nm  aspect {d['affine']['aspect_min_over_max']:.3f}  "
              f"rot {d['affine']['rotation_deg_(x_right_y_down_clockwise_positive)']:6.1f}  "
              f"mirror {d['affine']['mirrored']}  res {d['affine']['residual_nm_wrms']:5.1f} nm | "
              f"sim scale {d['similarity']['scale_nm_per_pitch']:5.1f} res {d['similarity']['residual_nm_wrms']:5.1f} | "
              f"perm-min {perm.min():5.1f}")
        b = d["huygens_conventions_best3"][0]
        print(f"{'':52s} Huygens best convention: k={b['factor_k']:.3f} res {b['residual_nm_wrms']:.1f} nm "
              f"(sign {b['rotation_sign']}, aspect on {b['aspect_on']}, yflip {b['y_flip']}, "
              f"rows vertical {b['layout_rows_vertical']}); centroid {np.round(d['intensity_centroid_pitch_units_xy'], 3)}")
    with open(args.json, "w") as f:
        json.dump(out, f, indent=1)
    print("wrote", args.json)
    plot(shifts_all, P, args.fig)
    print("wrote", args.fig)


def plot(shifts_all, P, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ink, ink2, surf = "#0b0b0b", "#52514e", "#fcfcfb"
    c_meas, c_model = "#2a78d6", "#eb6834"
    keys = [k for k in shifts_all if "Matrix 640 RAW" not in k]   # keep the figure readable
    n = len(keys)
    cols = 3
    rows = int(np.ceil(n / cols))
    fig, axs = plt.subplots(rows, cols, figsize=(4.2 * cols, 4.2 * rows), facecolor=surf)
    for ax, k in zip(np.ravel(axs), keys):
        S, w, A, t = shifts_all[k]
        M = P @ A.T + t
        ax.set_facecolor(surf)
        for e in range(N_EL):                  # measured -> fitted position of the same element
            ax.plot([S[e, 0], M[e, 0]], [S[e, 1], M[e, 1]], color="#a3a29c", lw=0.8, zorder=1)
        ax.scatter(M[:, 0], M[:, 1], s=90, facecolors="none", edgecolors=c_model, linewidths=2,
                   label="fitted hex lattice")
        ax.scatter(S[:, 0], S[:, 1], s=np.clip(w / w.max() * 120, 10, 120), color=c_meas,
                   edgecolors=surf, linewidths=1, label="measured shift", zorder=3)
        for e in range(N_EL):
            ax.annotate(str(e), S[e], xytext=(4, 4), textcoords="offset points", fontsize=7, color=ink2)
        ax.set_title(k, fontsize=9, color=ink)
        ax.set_aspect("equal")
        ax.invert_yaxis()                      # image convention: y (line index) down
        ax.axhline(0, color="#d8d7d2", lw=0.8, zorder=0)
        ax.axvline(0, color="#d8d7d2", lw=0.8, zorder=0)
        ax.tick_params(colors=ink2, labelsize=7)
        for sp in ax.spines.values():
            sp.set_color("#c3c2b7")
        ax.set_xlabel("x shift (nm)", fontsize=8, color=ink2)
        ax.set_ylabel("y shift (nm, down)", fontsize=8, color=ink2)
    for ax in np.ravel(axs)[n:]:
        ax.axis("off")
    h, l = np.ravel(axs)[0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower right", frameon=False, fontsize=9)
    fig.suptitle("Element image shifts vs sum of the other elements (labels = data element 0-22); "
                 "marker area ~ photons; grey line = residual to the fitted lattice", fontsize=10, color=ink)
    fig.tight_layout(rect=(0, 0.02, 1, 0.97))
    fig.savefig(path, dpi=110, facecolor=surf)


if __name__ == "__main__":
    main()
