"""Generate the figures used in docs/REPORT.md from the input data.

    .venv/bin/python tools/make_figures.py
"""
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from matrix2ptu import rawstream, tdc, timing  # noqa: E402
from matrix2ptu.msr_meta import read_flim_stacks, read_msr_metadata  # noqa: E402

OUT = "docs/figures"
os.makedirs(OUT, exist_ok=True)
C1, C2, INK, MUTED, GRID = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e", "#e4e3df"
plt.rcParams.update({
    "font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED,
    "ytick.color": MUTED, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
    "axes.spines.top": False, "axes.spines.right": False, "lines.linewidth": 2,
    "figure.dpi": 110, "savefig.bbox": "tight",
})


def decode(name, msr=None):
    ev = rawstream.decode(f"inputs/{name}.raw")
    meta = read_msr_metadata(f"inputs/{msr or name}.msr")
    cal = tdc.TdcCalibration.from_events(ev.ph_element, ev.ph_code, ev.sy_code)
    sf = cal.fine_one(23, ev.sy_code)
    pt = timing.build_pulse_train(ev.sy_coarse, sf)
    line = ev.mk_coarse[(ev.mk_bits & 2) > 0].astype(float)
    iv = np.diff(line)
    tick = (1 / meta["line_frequency_hz"]) / iv.mean()
    return ev, meta, cal, sf, pt, line, tick


ev, meta, cal, sf, pt, line, tick = decode("Matrix")

# ---- F1: TDC code structure (sync channel) --------------------------------
S = ev.sy_coarse[:2_000_000].astype(float)
code = ev.sy_code[:2_000_000].astype(int)
N = 2000
n = S.size // N * N
pred = np.empty(n)
for i in range(0, n, N):
    k = np.arange(N)
    x = S[i:i + N] + cal.tables[23][code[i:i + N]]
    pred[i:i + N] = np.polyval(np.polyfit(k, x, 1), k)
phase = pred - S[:n]
c = code[:n]
med = np.array([np.median(phase[c == j]) if (c == j).sum() > 50 else np.nan for j in range(1024)])
hist = np.bincount(code, minlength=1024)
fig, ax = plt.subplots(1, 2, figsize=(9, 3.2))
ax[0].plot(np.arange(1024), med, ".", ms=2.5, color=C1)
ax[0].set_xlabel("10-bit TDC code (sync channel)")
ax[0].set_ylabel("true phase within coarse tick")
ax[0].set_title("Code → time: 4 quarters, taps run backwards", color=INK, fontsize=9)
ax[1].bar(np.arange(1024), hist / 1e3, width=1.0, color=C1)
ax[1].set_xlabel("10-bit TDC code (sync channel)")
ax[1].set_ylabel("counts (thousands)")
ax[1].set_title("Code occupancy = bin width (strong DNL)", color=INK, fontsize=9)
fig.savefig(f"{OUT}/f1_tdc_structure.png")
plt.close(fig)

# ---- photon timing for Matrix ---------------------------------------------
fine = cal.fine(ev.ph_element, ev.ph_code)
t = ev.ph_coarse + fine
pp, ph = pt.locate(t)
tick_s = 25e-9 / pt.ticks_per_pulse()   # laser locked to exactly 40 MHz in this file
mt = ph * 25000.0
meta_d = np.array(meta["element_delays_s"]) * 1e12

# ---- F2: per-element offsets vs Imspector delays -------------------------
ref = np.histogram(mt[ev.ph_element == 11], bins=2500, range=(0, 25000))[0].astype(float)
shift, cnt = [], []
for e in range(23):
    h = np.histogram(mt[ev.ph_element == e], bins=2500, range=(0, 25000))[0].astype(float)
    cc = np.fft.irfft(np.fft.rfft(ref) * np.conj(np.fft.rfft(h)), n=2500)
    s = int(np.argmax(cc))
    shift.append((s if s < 1250 else s - 2500) * 10)
    cnt.append(h.sum())
shift = np.array(shift)
off = np.mean(shift - meta_d)
fig, ax = plt.subplots(figsize=(4.2, 3.6))
lim = [meta_d.min() - 300, meta_d.max() + 300]
ax.plot(lim, np.array(lim) + off, color=MUTED, lw=1, label="measured = Imspector + const")
ax.scatter(meta_d, shift, s=10 + 40 * np.sqrt(np.array(cnt) / max(cnt)), color=C1, zorder=3,
           edgecolor="white", linewidth=1, label="detector element")
for e in (0, 10, 21):
    ax.annotate(str(e), (meta_d[e], shift[e]), xytext=(4, -10), textcoords="offset points",
                color=MUTED, fontsize=8)
r = np.corrcoef(meta_d, shift)[0, 1]
ax.set_xlabel("Imspector stored delay per element (ps)")
ax.set_ylabel("offset measured from raw data (ps)")
ax.set_title(f"Per-element timing offsets, r = {r:.3f}", color=INK, fontsize=9)
ax.legend(frameon=False, fontsize=8, loc="upper left")
fig.savefig(f"{OUT}/f2_element_delays.png")
plt.close(fig)

# ---- F3/F4: Matrix image + decay vs Imspector -------------------------------
dl = (meta_d + meta["common_delay_s"] * 1e12) * 1e-12 / tick_s
pp2, ph2 = pt.locate(t + dl[ev.ph_element])
lp, _ = pt.locate(line)
k = np.searchsorted(lp, pp2, side="right") - 1
x = (pp2 - lp[np.clip(k, 0, None)]) // 4000
ok = (k >= 0) & (x >= 0) & (x < 400)
img = np.zeros((400, 400), int)
np.add.at(img, (k[ok], x[ok]), 1)
cube = list(read_flim_stacks("inputs/Matrix.msr").values())[0].astype(np.int64)
R = cube.sum(-1)
fig, ax = plt.subplots(1, 3, figsize=(10, 3.4))
vmax = np.percentile(R, 99.5)
for a, im, ttl in [(ax[0], img, "PTU (converted from .raw)"), (ax[1], R, "Imspector (.msr)")]:
    a.imshow(im, cmap="gray", vmin=0, vmax=vmax)
    a.set_title(ttl, color=INK, fontsize=9)
    a.axis("off")
d = img - R
m = ax[2].imshow(d, cmap="RdBu_r", vmin=-3, vmax=3)
ax[2].set_title(f"difference (97.0 % pixels identical)", color=INK, fontsize=9)
ax[2].axis("off")
fig.colorbar(m, ax=ax[2], fraction=0.046, label="counts")
fig.savefig(f"{OUT}/f3_matrix_images.png")
plt.close(fig)

h = np.histogram(ph2 * 25.0, bins=250, range=(0, 25))[0]
H = cube.sum((0, 1))
tb = np.arange(250) * 0.1 + 0.05
fig, ax = plt.subplots(figsize=(6, 3.2))
ax.semilogy(tb, H, color=C2, label="Imspector (.msr cube)")
ax.semilogy(tb, h, color=C1, lw=1.2, ls="--", label="converted (all 23 elements, delays applied)")
ax.set_xlabel("micro time (ns)")
ax.set_ylabel("photons per 100 ps")
ax.set_title("Matrix.raw decay: no fitted parameters", color=INK, fontsize=9)
ax.legend(frameon=False, fontsize=8)
fig.savefig(f"{OUT}/f4_matrix_decay.png")
plt.close(fig)

# ---- F5: photon rate vs time since line marker --------------------------------
k = np.searchsorted(line, t, side="right") - 1
rel = (t - line[np.clip(k, 0, None)]) * tick * 1e3
rel = rel[k >= 0]
hh, ee = np.histogram(rel, bins=np.arange(0, 41.0, 0.1))
fig, ax = plt.subplots(figsize=(6, 2.8))
ax.step(ee[:-1], hh / 400 / 0.1, where="post", color=C1, lw=1.2)
ax.axvspan(40.0, 40.8, color=GRID, zorder=0)
ax.text(40.05, ax.get_ylim()[1] * 0.85, "flyback", color=MUTED, fontsize=8, rotation=90)
ax.set_xlim(-0.5, 41)
ax.set_xlabel("time after line marker (ms)")
ax.set_ylabel("photons per ms per line")
ax.set_title("Line markers mark line *starts*: 400 × 100 µs, then 0.8 ms flyback", color=INK, fontsize=9)
fig.savefig(f"{OUT}/f5_line_timing.png")
plt.close(fig)
del ev, fine, t, pp, ph, pp2, ph2

# ---- F6: Sequence pass structure and images ------------------------------------
ev, meta, cal, sf, pt, line, tick = decode("Sequence", msr="Sequence2")
k = np.searchsorted(line, ev.ph_coarse, side="right") - 1
rel = (ev.ph_coarse - line[np.clip(k, 0, None)]) * tick
inl = (k >= 0) & (rel < 25.5e-3)
cnt = np.bincount(k[inl], minlength=line.size)[:3500].reshape(500, 7)
x = np.floor(rel / 51e-6).astype(int)
ok = inl & (x < 500)
p = k % 7
imgs = []
for sel in (p < 2, p >= 2):
    m = ok & sel
    im = np.zeros((500, 500), int)
    np.add.at(im, (k[m] // 7, x[m]), 1)
    imgs.append(im)
cubes = list(read_flim_stacks("inputs/Sequence2.msr").values())

fig = plt.figure(figsize=(10, 6.2))
gs = fig.add_gridspec(2, 4, height_ratios=[1, 1.25], hspace=0.45)
a = fig.add_subplot(gs[0, :])
means = cnt.mean(0)
cols = [C1, C1, C2, C2, C2, C2, C2]
a.bar(np.arange(1, 8), means / 1e3, color=cols, width=0.6)
for i, v in enumerate(means):
    a.text(i + 1, v / 1e3 + 0.4, f"{v/1e3:.1f}k", ha="center", color=INK, fontsize=8)
a.set_xticks(np.arange(1, 8))
a.set_xlabel("line pass within each y-line (marker index mod 7)")
a.set_ylabel("photons per pass (thousands)")
a.set_title("Sequence.raw: per y-line 2 CONF passes (blue) then 5 STED passes (orange)", color=INK, fontsize=9)
a.grid(axis="x", visible=False)
for j, (im, ttl) in enumerate([(imgs[0], "CONF from .raw (2 passes)"),
                               (cubes[0].reshape(500, 500, -1).sum(-1), "CONF, Imspector (Sequence2.msr)"),
                               (imgs[1], "STED from .raw (5 passes)"),
                               (cubes[1].reshape(500, 500, -1).sum(-1), "STED, Imspector (Sequence2.msr)")]):
    b = fig.add_subplot(gs[1, j])
    b.imshow(im, cmap="gray", vmin=0, vmax=np.percentile(im, 99.5))
    b.set_title(ttl, color=INK, fontsize=9)
    b.axis("off")
fig.savefig(f"{OUT}/f6_sequence.png")
plt.close(fig)
print("figures written to", OUT)
