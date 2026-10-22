"""Does FLOPA extract the same FLIM information from the PTU as Imspector stored?

Runs napari-flopa's own core functions (read_ptu_file + reconstruct_ptu_to_dataset,
no GUI) on a *pooled* PTU and compares, per pixel and per sequence, FLOPA's
photon count, mean arrival time and phasor (g, s) with the same quantities
computed from Imspector's FLIM cube in the .msr of the same acquisition
(100 ps bins, bin centres).  Pixels with fewer than --min-photons are skipped
for the lifetime quantities.

Usage:
    .venv/bin/python tools/flopa_flim_check.py outputs/MatrixOBF/MatrixOBF_pooled.ptu \
        inputs2/MatrixOBF.msr 125 125 1
    .venv/bin/python tools/flopa_flim_check.py outputs/OrangebeadsOBF/Sequence_orangebeadsOBF_pooled.ptu \
        inputs2/Sequence_orangebeadsOBF.msr 487 307 1,3
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from matrix2ptu.msr_meta import read_flim_stacks   # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ptu")
    ap.add_argument("msr")
    ap.add_argument("lines", type=int)
    ap.add_argument("pixels", type=int)
    ap.add_argument("acc")
    ap.add_argument("--min-photons", type=int, default=100)
    a = ap.parse_args()
    from napari_flopa.core.io.loader import read_ptu_file
    from napari_flopa.core.processing.reconstruction import reconstruct_ptu_to_dataset
    from tttrkit.ptuio.reconstructor import ScanConfig

    acc = tuple(int(x) for x in a.acc.split(","))
    data = read_ptu_file(a.ptu, header=False)
    cfg = ScanConfig(lines=a.lines, pixels=a.pixels, frames=1, line_accumulations=acc)
    ds = reconstruct_ptu_to_dataset(data, cfg, outputs=["photon_count", "mean_arrival_time", "phasor"])
    res_s = float(ds["tcspc_resolution"])
    omega = float(ds["omega"])                 # rad per TCSPC bin
    period_ns = 2 * np.pi / omega * res_s * 1e9
    cubes = read_flim_stacks(a.msr)
    names = list(cubes)
    print(f"FLOPA: tcspc resolution {res_s*1e12:g} ps, phasor period {period_ns:.4f} ns")
    ok_all = True
    for s in range(len(acc)):
        cnt = ds["photon_count"].values[0, s, :, :, 0].astype(float)
        mat = ds["mean_arrival_time"].values[0, s, :, :, 0].astype(float)
        g = ds["phasor_g"].values[0, s, :, :, 0].astype(float)
        sp = ds["phasor_s"].values[0, s, :, :, 0].astype(float)
        cube = cubes[names[s]].astype(float)            # (Y, X, T), 100 ps bins
        nT = cube.shape[-1]
        t = (np.arange(nT) + 0.5) * 0.1                 # ns
        n_ref = cube.sum(-1)
        with np.errstate(invalid="ignore", divide="ignore"):
            mat_ref = (cube * t).sum(-1) / n_ref
            w = 2 * np.pi / period_ns
            g_ref = (cube * np.cos(w * t)).sum(-1) / n_ref
            s_ref = (cube * np.sin(w * t)).sum(-1) / n_ref
        # FLOPA's mean arrival time unit: bins or ns -> decide from the data
        m = cnt >= a.min_photons
        scale = 1.0
        if np.nanmedian(mat[m]) > 3 * np.nanmedian(mat_ref[m]):
            scale = res_s * 1e9                          # bins -> ns
        mat_ns = mat * scale
        d_mat = (mat_ns - mat_ref)[m]
        d_g, d_s = (g - g_ref)[m], (sp - s_ref)[m]
        r_cnt = float(np.corrcoef(cnt.ravel(), n_ref.ravel())[0, 1])
        r_mat = float(np.corrcoef(mat_ns[m], mat_ref[m])[0, 1])
        print(f"sequence {s} vs Imspector '{names[s]}': {int(m.sum())} pixels with >= {a.min_photons} photons")
        print(f"  photon count: FLOPA {int(cnt.sum())} vs Imspector {int(n_ref.sum())}, r = {r_cnt:.5f}")
        print(f"  mean arrival time [ns] (FLOPA unit scale {scale:g}): median diff {np.median(d_mat)*1e3:+.1f} ps, "
              f"rms diff {np.sqrt(np.mean(d_mat**2))*1e3:.1f} ps, r = {r_mat:.4f}; "
              f"image mean {np.average(mat_ns[m], weights=cnt[m]):.4f} vs {np.average(mat_ref[m], weights=n_ref[m]):.4f}")
        print(f"  phasor g: median diff {np.median(d_g):+.4f}, rms {np.sqrt(np.mean(d_g**2)):.4f}; "
              f"s: median diff {np.median(d_s):+.4f}, rms {np.sqrt(np.mean(d_s**2)):.4f}")
        ok = r_cnt > 0.99 and abs(np.median(d_mat)) < 0.02 and abs(np.median(d_g)) < 0.01 \
            and abs(np.median(d_s)) < 0.01
        print(f"  -> {'PASS' if ok else 'FAIL'} (count r > 0.99, |median diff| < 20 ps and < 0.01 in g, s)")
        ok_all &= ok
    return 0 if ok_all else 1


if __name__ == "__main__":
    sys.exit(main())
