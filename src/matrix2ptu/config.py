"""Conversion settings: loading, defaults, and cross-checks.

Settings come from a YAML file (see configs/*.yaml).  Priority, as agreed with
the user:
  1. Values written in the YAML file are used as given (primary reference).
  2. Values set to "auto" (or omitted) are taken from the .msr metadata or
     measured from the .raw stream.
  3. Whatever the source, every value that can be checked against the .msr
     metadata or the raw stream is checked, and disagreements are reported as
     warnings (or errors when conversion would be wrong).

Keep these concepts separate (they are separate keys on purpose):
  * detector elements  (23 SPAD pixels; PTU channel field in resolved mode)
  * acquisition sequences / Imspector channels (e.g. CONF and STED)
  * line accumulations (repeated passes of the same y-line per sequence)
  * frames, lines (y), pixels (x)
"""
from __future__ import annotations

import copy
from dataclasses import dataclass, field

import yaml

AUTO = "auto"

DEFAULTS = {
    "input": {"raw": None, "msr": None},
    "output": {
        "directory": "outputs",
        "basename": None,
        "modes": ["detector_resolved", "pooled"],
        "split_sequences": False,
        "overwrite": False,
        # PTU HW_Type string. Honest default; if SymPhoTime refuses the file,
        # try the string Imspector writes for MultiHarp data ("MultiHarp 150 4P").
        "hw_type": "Abberior MATRIX (SPAD23) converted",
    },
    "scan": {
        "pixels_x": AUTO,
        "pixels_y": AUTO,
        "pixel_size_um": AUTO,
        "pixel_dwell_us": AUTO,
        # The dwell time applies to ONE pass of a line. With line
        # accumulation N, a pixel is exposed N x dwell in total for that
        # sequence (Imspector reports that product as the channel "dwelltime").
        "dwell_convention": "per_line_pass",
        "bidirectional": AUTO,
        "scan_type": "linear",
        "frames": 1,
        "line_start_offset_us": 0.0,
    },
    "acquisition": {
        "mode": "line_interleaved",
        "sequences": AUTO,   # list of {name, line_accumulation}
    },
    "detector": {
        "n_elements": 23,
        "sync_element": 23,
        "sync_divider": 32,
    },
    "timing": {
        "laser_frequency_hz": AUTO,
        "microtime_resolution_ps": 5.0,
        "apply_element_delays": True,
        "element_delays_ps": AUTO,
        "common_delay_ps": AUTO,
        "sync_smoothing_window": 64,
        "tdc_calibration_file": None,
    },
}


def _merge(base, over):
    out = copy.deepcopy(base)
    for k, v in (over or {}).items():
        if k not in out:
            raise ValueError(f"unknown configuration key '{k}'")
        if isinstance(out[k], dict) and isinstance(v, dict):
            for kk in v:
                if kk not in out[k]:
                    raise ValueError(f"unknown configuration key '{k}.{kk}'")
            out[k].update(v)
        else:
            out[k] = v
    return out


def load_config(path: str) -> dict:
    with open(path) as f:
        user = yaml.safe_load(f) or {}
    cfg = _merge(DEFAULTS, user)
    if not cfg["input"]["raw"]:
        raise ValueError("input.raw must be given")
    return cfg


@dataclass
class Resolved:
    """Final parameter set used for one conversion, with provenance."""
    values: dict = field(default_factory=dict)
    source: dict = field(default_factory=dict)
    warnings: list = field(default_factory=list)

    def set(self, key, value, source):
        self.values[key] = value
        self.source[key] = source

    def warn(self, msg):
        self.warnings.append(msg)

    def __getitem__(self, k):
        return self.values[k]


def _close(a, b, rel=1e-3):
    return abs(a - b) <= rel * max(abs(a), abs(b))


def resolve(cfg: dict, meta: dict | None) -> Resolved:
    """Combine user config with .msr metadata. Raw-stream checks happen later
    (they need the decoded data) in converter.check_against_raw()."""
    r = Resolved()
    sc, acq, tm, det = cfg["scan"], cfg["acquisition"], cfg["timing"], cfg["detector"]

    def pick(key, user_val, meta_val, unit_note="", rel=1e-3):
        if user_val != AUTO and user_val is not None:
            r.set(key, user_val, "config")
            if meta_val is not None:
                ok = (_close(float(user_val), float(meta_val), rel)
                      if isinstance(meta_val, (int, float)) else user_val == meta_val)
                if not ok:
                    r.warn(f"{key}: config value {user_val}{unit_note} differs from .msr "
                           f"metadata {meta_val}{unit_note}; using the config value")
        elif meta_val is not None:
            r.set(key, meta_val, ".msr metadata")
        else:
            raise ValueError(f"{key} is 'auto' but no .msr metadata is available; "
                             f"set scan/timing parameter '{key}' in the config file")

    m = meta or {}
    pick("pixels_x", sc["pixels_x"], m.get("pixels_x"))
    pick("pixels_y", sc["pixels_y"], m.get("pixels_y"))
    pick("pixel_size_um", sc["pixel_size_um"],
         m["pixel_size_x_m"] * 1e6 if "pixel_size_x_m" in m else None, " um")
    pick("pixel_dwell_us", sc["pixel_dwell_us"],
         m["pixel_dwell_per_pass_s"] * 1e6 if "pixel_dwell_per_pass_s" in m else None, " us")
    pick("bidirectional", sc["bidirectional"], m.get("bidirectional"))
    r.set("frames", int(sc["frames"]), "config")
    r.set("scan_type", sc["scan_type"], "config")
    r.set("line_start_offset_us", float(sc["line_start_offset_us"]), "config")
    if sc["dwell_convention"] != "per_line_pass":
        raise ValueError("only dwell_convention: per_line_pass is implemented")
    r.set("dwell_convention", "per_line_pass", "config")

    if r["bidirectional"]:
        raise ValueError("bidirectional scanning is not implemented (no test data); "
                         "the MATRIX stream gives only line-start markers")
    if r["scan_type"] != "linear":
        raise ValueError("only linear (non-harmonic) scanning is implemented")

    # sequences
    msr_seq = None
    if meta and meta.get("channels"):
        steps = meta["line_steps"]
        chans = meta["channels"]
        # Assign each Imspector channel to the line step(s) it is active in.
        # Only the simple case "one channel per step" is mapped automatically.
        if len(chans) == len(steps):
            msr_seq = [{"name": c["name"], "line_accumulation": s} for c, s in zip(chans, steps)]
            for c, s in zip(chans, steps):
                rep = c["dwelltime_reported_s"]
                per = m["pixel_dwell_per_pass_s"]
                if per and rep and not _close(rep, per * s, 1e-6):
                    r.warn(f"sequence '{c['name']}': Imspector reports dwell {rep*1e6:g} us, "
                           f"expected per-pass dwell x accumulation = {per*s*1e6:g} us")
        else:
            r.warn(f".msr has {len(chans)} channels but {len(steps)} line steps; "
                   "sequence mapping must be given in the config")
    if acq["sequences"] in (AUTO, None):
        if msr_seq is None:
            raise ValueError("acquisition.sequences is 'auto' but cannot be derived from the "
                             ".msr file; list the sequences with their line_accumulation")
        r.set("sequences", msr_seq, ".msr metadata")
    else:
        seqs = [{"name": str(s.get("name", f"seq{i}")),
                 "line_accumulation": int(s["line_accumulation"])}
                for i, s in enumerate(acq["sequences"])]
        r.set("sequences", seqs, "config")
        if msr_seq is not None:
            a = [s["line_accumulation"] for s in seqs]
            b = [s["line_accumulation"] for s in msr_seq]
            if a != b:
                r.warn(f"line accumulations in config {a} differ from .msr line steps {b}")
    if any(s["line_accumulation"] < 1 for s in r["sequences"]):
        raise ValueError("line_accumulation must be >= 1")
    if acq["mode"] != "line_interleaved":
        raise ValueError("acquisition.mode must be 'line_interleaved' (a single sequence is "
                         "the special case with one entry)")
    r.set("acquisition_mode", "line_interleaved", "config")

    # detector
    r.set("n_elements", int(det["n_elements"]), "config")
    r.set("sync_element", int(det["sync_element"]), "config")
    r.set("sync_divider", int(det["sync_divider"]), "config")
    if meta and meta.get("array_id") and meta["array_id"] != f"hex_{det['n_elements']}":
        r.warn(f".msr array_id is {meta['array_id']}, config n_elements is {det['n_elements']}")

    # timing
    r.set("microtime_resolution_ps", float(tm["microtime_resolution_ps"]), "config")
    r.set("sync_smoothing_window", int(tm["sync_smoothing_window"]), "config")
    r.set("tdc_calibration_file", tm["tdc_calibration_file"], "config")
    r.set("apply_element_delays", bool(tm["apply_element_delays"]), "config")
    if tm["apply_element_delays"]:
        ed = tm["element_delays_ps"]
        if ed in (AUTO, None):
            if not m.get("element_delays_s"):
                raise ValueError("element_delays_ps is 'auto' but no .msr metadata; give 23 "
                                 "values in ps, or set apply_element_delays: false")
            r.set("element_delays_ps", [d * 1e12 for d in m["element_delays_s"]], ".msr metadata")
        else:
            r.set("element_delays_ps", [float(v) for v in ed], "config")
        if len(r["element_delays_ps"]) != r["n_elements"]:
            raise ValueError(f"{len(r['element_delays_ps'])} element delays given for "
                             f"{r['n_elements']} elements")
        cd = tm["common_delay_ps"]
        if cd in (AUTO, None):
            if "common_delay_s" not in m:
                raise ValueError("common_delay_ps is 'auto' but no .msr metadata")
            r.set("common_delay_ps", m["common_delay_s"] * 1e12, ".msr metadata")
        else:
            r.set("common_delay_ps", float(cd), "config")
    else:
        r.set("element_delays_ps", [0.0] * r["n_elements"], "disabled")
        cd = tm["common_delay_ps"]
        r.set("common_delay_ps", 0.0 if cd in (AUTO, None) else float(cd),
              "disabled" if cd in (AUTO, None) else "config")

    # laser frequency: resolved later against the raw stream; keep user value
    r.set("laser_frequency_hz", tm["laser_frequency_hz"], "config" if tm["laser_frequency_hz"] != AUTO else "pending")
    r.set("line_frequency_hz_msr", m.get("line_frequency_hz"), ".msr metadata" if m else "n/a")

    # geometry for header
    r.set("x0_um", (m["offset_x_m"] - m["fov_x_m"] / 2) * 1e6 if m else 0.0, ".msr metadata" if m else "default")
    r.set("y0_um", (m["offset_y_m"] - m["fov_y_m"] / 2) * 1e6 if m else 0.0, ".msr metadata" if m else "default")
    if m.get("crosstalk_filter_active"):
        r.warn("Imspector crosstalk filter was active; the raw stream is unfiltered")
    for c in (m.get("channels") or []):
        if c.get("gating_active"):
            r.warn(f"Imspector time gating was active for '{c['name']}' "
                   f"(delay {c['gating_delay_s']*1e9:g} ns, width {c['gating_width_s']*1e9:g} ns). "
                   "The PTU contains all photons (no gating applied); apply the gate in analysis if wanted.")
    return r
