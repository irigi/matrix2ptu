"""Headless check using napari-flopa's own core I/O + reconstruction functions
(no napari GUI). Usage:
    .venv/bin/python tools/flopa_headless_check.py outputs/Sequence/Sequence_detector_resolved.ptu 500 500 2,5
"""
import sys

import numpy as np
from napari_flopa.core.io.loader import read_ptu_file
from napari_flopa.core.processing.reconstruction import reconstruct_ptu_to_dataset
from tttrkit.ptuio.reconstructor import ScanConfig

path, lines, pixels, acc = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
acc = tuple(int(a) for a in acc.split(","))
data = read_ptu_file(path, header=False)
cfg = ScanConfig(lines=lines, pixels=pixels, frames=1, line_accumulations=acc)
ds = reconstruct_ptu_to_dataset(data, cfg, outputs=["photon_count", "mean_arrival_time", "phasor"])
pc = ds["photon_count"].values
print("dataset dims:", dict(ds.sizes))
for s in range(pc.shape[1]):
    img = pc[0, s].sum(-1)
    print(f"sequence {s}: photons {int(img.sum())}, active channels {int((pc[0, s].sum((0, 1)) > 0).sum())}")
g = ds["phasor_g"].values
print("phasor g finite fraction:", float(np.isfinite(g).mean()))
