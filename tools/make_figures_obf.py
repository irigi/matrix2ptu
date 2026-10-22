"""Figures for docs/OBF_REPORT.md, made from the converted PTU files.

    .venv/bin/python tools/make_figures_obf.py
(run the four conversions in configs/{matrix_obf,matrix_raw640,orangebeads_obf,orangebeads_raw}.yaml first)
"""
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = os.path.dirname(__file__)
sys.path.insert(0, os.path.join(HERE, "..", "src"))
sys.path.insert(0, HERE)
from compare_pair import load, register  # noqa: E402
from matrix2ptu.msr_meta import read_flim_stacks  # noqa: E402
from matrix2ptu.validate import iter_ptu_records  # noqa: E402

OUT = os.path.join(HERE, "..", "docs", "figures")
C1, C2, C3, INK, MUTED, GRID = "#2a78d6", "#eb6834", "#0b0b0b", "#0b0b0b", "#52514e", "#e4e3df"
plt.rcParams.update({
    "font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED,
    "ytick.color": MUTED, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
    "axes.spines.top": False, "axes.spines.right": False, "lines.linewidth": 1.6,
    "figure.dpi": 110, "savefig.bbox": "tight",
})
O = os.path.join(HERE, "..", "outputs")
I2 = os.path.join(HERE, "..", "inputs2")


def to100(h, res_ps=5):
    f = int(round(100 / res_ps))
    n = h.size // f * f
    return h[:n].reshape(-1, f).sum(1)


def img_panel(ax, im, title):
    ax.imshow(im, cmap="gray", vmax=np.percentile(im, 99.5), interpolation="nearest")
    ax.set_title(title, fontsize=8)
    ax.set_xticks([]); ax.set_yticks([]); ax.grid(False)


def fig_pairs():
    pairs = [("Dataset A (640 nm)", "MatrixRAW/MatrixRAW", "MatrixOBF/MatrixOBF",
              "MatrixRAW.msr", "MatrixOBF.msr"),
             ("Dataset B (561 nm)", "OrangebeadsRAW/Sequence_orangebeadsRAW",
              "OrangebeadsOBF/Sequence_orangebeadsOBF",
              "Sequence_orangebeadsRAW.msr", "Sequence_orangebeadsOBF.msr")]
    fig = plt.figure(figsize=(11, 8.2))
    gs = fig.add_gridspec(3, 4)
    row = 0
    for title, pa, pb, ma, mb in pairs:
        A, B = load(os.path.join(O, pa)), load(os.path.join(O, pb))
        ca, cb = read_flim_stacks(os.path.join(I2, ma)), read_flim_stacks(os.path.join(I2, mb))
        for i, s in enumerate(A["seqs"]):
            r = row + i
            img_panel(fig.add_subplot(gs[r, 0]), A["img"][i], f"{title} {s['name']}\nRAW acquisition (PTU)")
            dy, dx, rr = register(A["img"][i].astype(float), B["img"][i].astype(float))
            img_panel(fig.add_subplot(gs[r, 1]), B["img"][i],
                      f"OBF acquisition (PTU)\nshift vs RAW ({dy}, {dx}) px, r = {rr:.3f}")
            ax = fig.add_subplot(gs[r, 2:])
            t = (np.arange(250) + 0.5) * 0.1
            for h, lab, col, ls in ((to100(A["dec"][i]), "RAW → PTU", C1, "-"),
                                    (to100(B["dec"][i]), "OBF → PTU", C2, "-"),
                                    (list(cb.values())[i].reshape(-1, 250).sum(0), "Imspector cube (OBF .msr)", C3, ":")):
                h = np.asarray(h, float)[:250]
                ax.semilogy(t[:h.size], h / h.sum(), ls, color=col, label=lab, lw=1.4 if ls == "-" else 1.2)
            ax.set_xlim(0, 25)
            ax.set_ylim(1e-5, None)
            ax.set_xlabel("micro time (ns)")
            ax.set_ylabel("normalised counts / 100 ps")
            ax.set_title(f"{s['name']}: decays (image photons)", fontsize=8)
            if r == 0:
                ax.legend(frameon=False, fontsize=8)
        row += len(A["seqs"])
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "f8_raw_vs_obf.png"))
    plt.close(fig)


def fig_pass_timing():
    """Photon rate vs time after each line marker, per pass position (dataset B)."""
    path = os.path.join(O, "OrangebeadsOBF/Sequence_orangebeadsOBF_pooled.ptu")
    ls = np.concatenate([t[(d["marker"] & 1) > 0] for _, d, t in iter_ptu_records(path)])
    le = np.concatenate([t[(d["marker"] & 2) > 0] for _, d, t in iter_ptu_records(path)])
    period = int(np.median(np.diff(ls)))
    L = int(np.median(le - ls))
    nb = 200
    H = np.zeros((4, nb))
    for _, d, t in iter_ptu_records(path):
        tp = t[(d["channel"] >= 0) & (d["marker"] == 0)]
        k = np.searchsorted(ls, tp, side="right") - 1
        ok = k >= 0
        rel = tp[ok] - ls[k[ok]]
        pos = k[ok] % 4
        m = rel < period
        np.add.at(H, (pos[m], (rel[m] * nb // period).astype(np.int64)), 1)
    Tp_ms = 1 / 40000568 * 1e3
    tt = (np.arange(nb) + 0.5) * period / nb * Tp_ms
    fig, ax = plt.subplots(figsize=(8, 3.2))
    labels = ["pass 1: CONF (EXC561)", "pass 2: STED", "pass 3: STED", "pass 4: STED"]
    cols = [C1, C2, "#1baf7a", "#6250d6"]
    n_lines = ls.size // 4
    for i in range(4):
        ax.plot(tt, H[i] / n_lines / (period / nb * Tp_ms * 1e-3) / 1e6, color=cols[i], label=labels[i])
    ax.axvspan(L * Tp_ms, period * Tp_ms, color=GRID, alpha=0.8, lw=0)
    ax.text(L * Tp_ms - 0.1, 0.02, "flyback (outside image) →", fontsize=8, color=MUTED, ha="right", va="bottom")
    ax.set_xlabel("time after line marker (ms)")
    ax.set_ylabel("photon rate (MHz), mean over y-lines")
    ax.set_title("Dataset B (OBF): photon rate within one line pass, per pass of the 4-pass cycle\n(the three STED passes overlap)", fontsize=9)
    ax.legend(frameon=False, fontsize=8, ncol=2, loc="upper left")
    fig.savefig(os.path.join(OUT, "f9_pass_timing.png"))
    plt.close(fig)


if __name__ == "__main__":
    if "--only-timing" not in sys.argv:
        fig_pairs()
    fig_pass_timing()
    print("wrote f8_raw_vs_obf.png, f9_pass_timing.png")
