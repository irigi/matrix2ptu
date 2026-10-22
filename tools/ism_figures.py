"""Figures for docs/ISM_REPORT.md (run tools/ism_inputs3.py stages first)."""
import json
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import tifffile  # noqa: E402
from matplotlib.patches import Ellipse  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
R = os.path.join(ROOT, "outputs", "inputs3_reassignment")
OUT = os.path.join(ROOT, "docs", "figures")
C1, C2, C3, C4, INK, MUTED, GRID = "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#0b0b0b", "#52514e", "#e4e3df"
plt.rcParams.update({
    "font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED,
    "ytick.color": MUTED, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
    "axes.spines.top": False, "axes.spines.right": False, "lines.linewidth": 1.6,
    "figure.dpi": 110, "savefig.bbox": "tight",
})


def load(p):
    with open(os.path.join(R, p)) as f:
        return json.load(f)


def img(ax, im, title, vmax=None):
    ax.imshow(im, cmap="magma", vmin=0, vmax=vmax if vmax is not None else np.percentile(im, 99.7),
              interpolation="nearest")
    ax.set_title(title, fontsize=8)
    ax.set_xticks([]); ax.set_yticks([]); ax.grid(False)


def fig_beads(cal):
    fig, axs = plt.subplots(1, 2, figsize=(10, 5.2))
    for ax, s in zip(axs, ("CONF", "STED")):
        im = tifffile.imread(os.path.join(R, "3_reassigned", "Beads", "images",
                                          f"Beads_seq{0 if s == 'CONF' else 1}-{s}_pooled_original.tif"))
        img(ax, im, f"Beads {s}, pooled, before reassignment")
        c = cal["sequences"][s]
        used = {(b["x"], b["y"]) for b in c["beads_used"]}
        crowd = {(b["x"], b["y"]) for b in c["beads_excluded_crowded"]}
        out = {(b["x"], b["y"]) for b in c["beads_rejected_outlier"]}
        for b in c["bead_table"]:
            p = (b["x"], b["y"])
            if p in used:
                col, lab = C3, "used for calibration"
            elif p in crowd:
                col, lab = C4, "isolated in STED but >1 neighbour in fit window"
            elif p in out:
                col, lab = C2, "outlier (rejected)"
            else:
                col, lab = "#9b9a96", "not isolated / edge / not single in STED"
            ax.plot(b["x"], b["y"], "o", mfc="none", mec=col, ms=11, mew=1.4, label=lab)
        h, l = ax.get_legend_handles_labels()
        u = dict(zip(l, h))
        ax.legend(u.values(), u.keys(), fontsize=7, loc="lower center", bbox_to_anchor=(0.5, -0.2), frameon=False, ncol=2)
    fig.savefig(os.path.join(OUT, "f10_bead_selection.png"))
    plt.close(fig)


def fig_shifts(cal, tr):
    from matrix2ptu.beadcal import layout
    fig, axs = plt.subplots(1, 3, figsize=(13, 4.6), gridspec_kw={"width_ratios": [1, 1, 0.8]})
    for ax, s in zip(axs[:2], ("CONF", "STED")):
        c = cal["sequences"][s]
        S = np.asarray(c["shift_nm"])
        U = np.asarray(c["uncertainty_nm"])
        ax.axhline(0, color=MUTED, lw=0.6); ax.axvline(0, color=MUTED, lw=0.6)
        for e in range(23):
            ax.add_patch(Ellipse(S[e], 2 * 2 * U[e, 0], 2 * 2 * U[e, 1], fc=C1, alpha=0.18, ec="none"))
        ax.plot(S[:, 0], S[:, 1], "o", color=C1, ms=4, label="bead calibration (Beads), ±2σ")
        if s == "CONF":
            A = np.asarray(tr["CONF"]["Dense"]["self_estimate_shift_px"]) * 20
            ax.plot(A[:, 0], A[:, 1], "^", color=C2, ms=4, label="Dense self-estimate (cross-correlation)")
            for e in range(23):
                ax.plot([S[e, 0], A[e, 0]], [S[e, 1], A[e, 1]], color=C2, lw=0.5, alpha=0.6)
        for e in range(23):
            ax.annotate(str(e + 1), S[e], textcoords="offset points", xytext=(3, 3), fontsize=6.5, color=INK)
        ax.set_aspect("equal")
        ax.invert_yaxis()
        ax.set_xlabel("shift x (nm), + = later in the line")
        ax.set_ylabel("shift y (nm), + = later line (down)")
        ax.set_title(f"{s}: element image position relative to pooled image\n(labels = Imspector element 1–23)", fontsize=8)
        ax.legend(fontsize=7, frameon=False, loc="lower left")
    P = layout()
    ax = axs[2]
    ax.plot(P[:, 0], P[:, 1], "o", color=MUTED, ms=14, mfc="none")
    for e in range(23):
        ax.text(P[e, 0], P[e, 1], str(e + 1), ha="center", va="center", fontsize=7)
    ax.set_aspect("equal"); ax.invert_yaxis(); ax.grid(False)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_title("Imspector element layout (labels)\nas drawn on screen", fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "f11_element_shifts.png"))
    plt.close(fig)


def fig_loo(cal, val):
    fig, axs = plt.subplots(1, 3, figsize=(13, 3.8))
    ax = axs[0]
    for s, col in (("CONF", C1), ("STED", C2)):
        sm = cal["sequences"][s]["loo_validation"]["summary"]
        ks = sorted(float(k[1:]) for k in sm)
        m = [sm[f"k{k:g}"]["fwhm_ratio"]["mean"] for k in ks]
        e = [sm[f"k{k:g}"]["fwhm_ratio"]["sem"] for k in ks]
        ax.errorbar(ks, m, yerr=e, color=col, marker="o", ms=4, capsize=2, label=f"{s} (n = {sm['k1']['fwhm_ratio']['n']} beads)")
    ax.axhline(1, color=MUTED, lw=0.8, ls=":")
    ax.axvline(1, color=MUTED, lw=0.8, ls=":")
    ax.set_xlabel("shift scale k (1 = measured shifts)")
    ax.set_ylabel("FWHM after / before (mean ± s.e.m.)")
    ax.set_title("Held-out beads (leave-one-bead-out)", fontsize=8)
    ax.legend(frameon=False, fontsize=7)
    for ax, s in zip(axs[1:], ("CONF", "STED")):
        rows = val[f"Beads_{s}"]["object_fwhm"]
        for nm, col, mk in (("calibration_beads", C1, "o"), ("crowded_or_overlapping", C4, "s")):
            if nm not in rows:
                continue
            b = np.array([x["before"]["sigma_px"] for x in rows[nm]["per_object"]]) * 47.096
            a = np.array([x["after_events"]["sigma_px"] for x in rows[nm]["per_object"]]) * 47.096
            ax.plot(b, a, mk, color=col, ms=4, mfc="none" if mk == "s" else col,
                    label=f"{nm.replace('_', ' ')} (n = {len(b)})")
        lo, hi = ax.get_xlim()
        ax.plot([lo, hi], [lo, hi], color=MUTED, lw=0.8, ls=":")
        ax.set_xlabel("FWHM before (nm)"); ax.set_ylabel("FWHM after, event-level PTU (nm)")
        ax.set_title(f"Beads {s}: per-bead FWHM (Gaussian fit)", fontsize=8)
        ax.legend(frameon=False, fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "f12_bead_validation.png"))
    plt.close(fig)


def fig_before_after():
    fig, axs = plt.subplots(4, 4, figsize=(12, 12.4))
    crop = (slice(20, 140), slice(20, 140))
    for r, (D, s, i) in enumerate((("Beads", "CONF", 0), ("Beads", "STED", 1), ("Dense", "CONF", 0), ("Dense", "STED", 1))):
        base = os.path.join(R, "3_reassigned", D, "images", f"{D}_seq{i}-{s}")
        b = tifffile.imread(base + "_pooled_original.tif")[crop]
        a = tifffile.imread(base + "_pooled_reassigned-events.tif")[crop]
        el = tifffile.imread(base + "_elements_reassigned-events.tif")
        vmax = np.percentile(np.r_[b.ravel(), a.ravel()], 99.8)
        img(axs[r, 0], b, f"{D} {s}: before (pooled PTU)", vmax)
        img(axs[r, 1], a, f"{D} {s}: after (reassigned pooled PTU)", vmax)
        d = a - b
        lim = np.percentile(np.abs(d), 99.5)
        axs[r, 2].imshow(d, cmap="RdBu_r", vmin=-lim, vmax=lim, interpolation="nearest")
        axs[r, 2].set_title("after − before (red = gained)", fontsize=8)
        axs[r, 2].set_xticks([]); axs[r, 2].set_yticks([]); axs[r, 2].grid(False)
        # line profile through the brightest object in the crop
        y, x = np.unravel_index(np.argmax(b), b.shape)
        ax = axs[r, 3]
        xs = (np.arange(b.shape[1]) - x) * 20
        ax.plot(xs, b[y], color=MUTED, label="before")
        ax.plot(xs, a[y], color=C1, label="after")
        ax.set_xlim(-600 if s == "CONF" else -300, 600 if s == "CONF" else 300)
        ax.set_xlabel("x (nm) through the brightest pixel row")
        ax.set_ylabel("counts / pixel")
        ax.legend(frameon=False, fontsize=7)
        ax.set_title(f"{D} {s}: line profile", fontsize=8)
        axs[r, 0].axhline(y, color="w", lw=0.4, alpha=0.6)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "f13_before_after.png"))
    plt.close(fig)


def fig_transfer(tr):
    fig, axs = plt.subplots(1, 2, figsize=(11, 3.8))
    for ax, s in zip(axs, ("CONF", "STED")):
        for D, ls in (("Beads", "-"), ("Dense", "--")):
            bg = tr[s][D]["band_gain"]
            per = [np.mean(b["band_period_nm"]) for b in bg["bands"]]
            sig = [b["significant"] for b in bg["bands"]]
            for lab, col in (("bead calibration", C1), ("1.4 x bead", C2), ("self-estimate", C3)):
                g = np.array([np.nan if v is None else v for v in bg["gains"][lab]])
                g[~np.array(sig)] = np.nan
                ax.plot(per, g, ls, marker="o", ms=3.5, color=col, label=f"{D}: {lab}")
        ax.axhline(1, color=MUTED, lw=0.8, ls=":")
        ax.set_xscale("log")
        ax.invert_xaxis()
        ax.set_xlabel("spatial period (nm), finer →")
        ax.set_ylabel("signal amplitude gain vs. no reassignment")
        ax.set_title(f"{s}: gain per frequency band (solid = Beads, dashed = Dense;\nbands with significant signal only)", fontsize=8)
        ax.legend(frameon=False, fontsize=6.5, ncol=2)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "f14_transfer_gain.png"))
    plt.close(fig)


def main():
    import sys
    sys.path.insert(0, os.path.join(ROOT, "src"))
    cal = load("2_calibration/Beads_ism_calibration.json")
    tr = load("4_validation/transfer_assessment.json")
    val = load("4_validation/validation.json")
    fig_beads(cal)
    fig_shifts(cal, tr)
    fig_loo(cal, val)
    fig_before_after()
    fig_transfer(tr)
    print("wrote f10-f14 in docs/figures")


if __name__ == "__main__":
    main()
