"""Extract acquisition parameters relevant to MATRIX .raw decoding from an
Imspector .msr/.obf file (same measurement).

The Imspector parameter tree is stored as XML in each stack's footer.  It is
flattened into "path = value" strings and a handful of known keys are read.
Everything returned here is *advisory*: the converter compares it with the
user configuration and with what the .raw stream itself shows, and reports any
disagreement.
"""
from __future__ import annotations

import xml.etree.ElementTree as ET

import numpy as np

from .obf import read_obf


def flatten_xml(text: str) -> dict[str, str]:
    root = ET.fromstring(text)
    out: dict[str, str] = {}

    def walk(e, path):
        kids = list(e)
        if not kids:
            out[path] = (e.text or "").strip()
        for i, k in enumerate(kids):
            name = f"item[{i}]" if k.tag == "item" else k.tag
            walk(k, f"{path}/{name}")

    walk(root, "")
    return out


def _floats(s):
    return [float(v) for v in s.split()] if s else []


def read_msr_metadata(path: str) -> dict:
    stacks, _file_docs, buf = read_obf(path)
    params = None
    for st in stacks:
        for doc in st.xml_docs:
            if "<ExpControl" in doc and "<measurement" in doc:
                params = flatten_xml(doc)
                break
        if params:
            break
    if params is None:
        raise ValueError(f"{path}: no Imspector measurement parameters found")

    g = params.get
    ms = "/doc/ExpControl/measurement"
    sc = "/doc/ExpControl/scan/range"

    # which Matrix stream was active/streaming
    streams = []
    for name in ("Matrix1", "Matrix2"):
        if g(f"/doc/{name}/data/streaming/enable") == "1":
            streams.append(name)
    det = streams[0] if streams else "Matrix1"

    channels = []
    i = 0
    while g(f"{ms}/channels/item[{i}]/name") is not None:
        p = f"{ms}/channels/item[{i}]"
        lasers = []
        j = 0
        while g(f"{p}/lasers/item[{j}]/id") is not None:
            if g(f"{p}/lasers/item[{j}]/active") == "1":
                lasers.append((g(f"{p}/lasers/item[{j}]/id"),
                               float(g(f"{p}/lasers/item[{j}]/power/calibrated") or "nan")))
            j += 1
        linestep_active = [int(v) for v in (g(f"{p}/steps/linestep_active") or "").split()]
        channels.append({
            "name": g(f"{p}/name"),
            "detector": g(f"{p}/detsel/detector"),
            "stream": g(f"{p}/tcspc/stream"),
            "dwelltime_reported_s": float(g(f"{p}/dwelltime") or "nan"),
            "linestep_active": linestep_active,
            "gating_active": g(f"{p}/gating/active") == "1",
            "gating_delay_s": float(g(f"{p}/gating/delay") or "nan"),
            "gating_width_s": float(g(f"{p}/gating/width") or "nan"),
            "flim_bin_s": float(g(f"{p}/tcspc/flim_psz_used") or "nan"),
            "flim_bins": int(float(g(f"{p}/tcspc/flim_res_used") or "0")),
            "element_mask": g(f"{p}/tcspc/mc_mask"),
            "lasers_active": lasers,
        })
        i += 1

    ls_dur = [int(float(v)) for v in (g(f"{ms}/line_steps/step_duration") or "").split()]
    ls_act = [int(v) for v in (g(f"{ms}/line_steps/step_active") or "").split()]
    line_steps = [d for d, a in zip(ls_dur, ls_act) if a] if g(f"{ms}/line_steps/active") == "1" else [1]

    meta = {
        "source": path,
        "imspector_version": g("/spec_version") or g("/item[1]/spec_version"),
        "detector_stream": det,
        "firmware": g(f"/hwr/{det}/firmware/bitfile_id"),
        "array_id": g(f"/hwr/{det}/firmware/array_id"),
        "element_delays_s": _floats(g(f"/hwr/{det}/delays/channel")),
        "common_delay_s": float(g(f"/hwr/{det}/delays/common") or "0"),
        "crosstalk_filter_active": g(f"/doc/{det}/crosstalk_filter/active") == "1",
        "sync_frequency_hz": float(g("/doc/ExpControl/timing/sync_source/frequency") or "nan"),
        "sync_source_input": g("/doc/ExpControl/timing/sync_source/input"),
        # Imspector's measured frequencies of the (external) sync inputs
        "ext_frequencies_hz": _floats(g("/doc/ExpControl/timing/ext_frequencies")),
        "pixels_x": int(float(g(f"{sc}/x/res") or 0)),
        "pixels_y": int(float(g(f"{sc}/y/res") or 0)),
        "pixel_size_x_m": float(g(f"{sc}/x/psz") or "nan"),
        "pixel_size_y_m": float(g(f"{sc}/y/psz") or "nan"),
        "fov_x_m": float(g(f"{sc}/x/len") or "nan"),
        "fov_y_m": float(g(f"{sc}/y/len") or "nan"),
        "offset_x_m": float(g(f"{sc}/x/off") or "nan"),
        "offset_y_m": float(g(f"{sc}/y/off") or "nan"),
        "scan_mode": g(f"{sc}/scanmode"),
        "bidirectional": g("/doc/Measurement/axes/bidirectional_scan") == "1",
        "line_frequency_hz": float(g("/doc/ExpControl/scan/line_frequency") or "nan"),
        "pixel_dwell_per_pass_s": float(g(f"{ms}/pixel_steps/total") or "nan"),
        "pixel_step_repetitions": int(float(g(f"{ms}/pixel_steps/repetitions") or 1)),
        "line_steps": line_steps,
        "line_step_repetitions": int(float(g(f"{ms}/line_steps/repetitions") or 1)),
        "channels": channels,
        "stacks": [{"name": s.name, "shape": s.shape} for s in stacks],
    }
    return meta


def read_flim_stacks(path: str) -> dict[str, np.ndarray]:
    """Return {stack name: array (Y, X, T)} for Imspector FLIM stacks (validation)."""
    stacks, _, buf = read_obf(path)
    out = {}
    for s in stacks:
        if len(s.shape) == 4 and s.shape[-1] > 1 and s.shape[1] > 1:
            out[s.name] = s.read_data(buf)[0]
    return out
