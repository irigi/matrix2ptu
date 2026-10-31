"""MATRIX .raw -> PTU conversion pipeline.

Steps (each documented in docs/REPORT.md):
  1. decode the 32-bit stream, unwrap the coarse counter      (rawstream)
  2. calibrate each TDC by code density                       (tdc)
  3. rebuild the laser pulse train from the /32 sync records  (timing)
  4. per photon: absolute time -> optional per-element delay correction
     -> laser pulse index (PTU nsync) + microtime (PTU dtime)
  5. markers: frame + line-start from the stream; line-stop synthesised at
     line start + pixels_x * dwell (the stream has no line-end marker, and
     PTU readers need one)
  6. write PTU (HydraHarp-v2 T3 records), detector-resolved and/or pooled,
     optionally one file per acquisition sequence.
"""
from __future__ import annotations

import datetime as dt
import json
import os
import time

import numpy as np

from . import __version__, rawstream, tdc, timing
from .config import AUTO, Resolved, load_config, resolve, source_path
from .msr_meta import compare_metadata, read_msr_metadata
from .ptu_writer import REC_TYPE_HYDRAHARP2_T3, WRAP, new_guid, write_ptu_t3

MARK_LINE_START = 1   # marker bit values (ImgHdr_LineStart = 1 -> bit value 1)
MARK_LINE_STOP = 2    # ImgHdr_LineStop = 2 -> bit value 2
MARK_FRAME = 4        # ImgHdr_Frame = 3 -> bit value 4


class ConversionError(RuntimeError):
    pass


def log(msg):
    print(f"[matrix2ptu] {msg}", flush=True)


# --------------------------------------------------------------------------- #
# Decoding and timing (shared by conversion and validation)
# --------------------------------------------------------------------------- #
class DecodedAcquisition:
    """All per-event quantities derived from one event stream (.raw or streamed .obf).

    Per photon (sorted by laser pulse; photons with the same pulse keep their
    order in the stream): `pulse` (int64, PTU nsync), `dtime` (uint16),
    `element` (uint8), `pass_of_photon` (int32, line pass, -1 before the first
    line marker).  With keep_validation_data=True also `t` (absolute time in
    coarse ticks, after delay correction) and `phase` (fraction of the laser
    period); these are only needed by validate.py and cost 16 bytes/photon.
    """

    BLOCK = 1_000_000     # photons per block (bounds temporary memory)

    def __init__(self, cfg: dict, r: Resolved, meta: dict | None,
                 keep_validation_data: bool = False):
        self.r = r
        raw = source_path(cfg)
        t0 = time.time()
        log(f"decoding {raw}")
        ev = rawstream.decode(raw, sync_element=r["sync_element"], keep_sync=False)
        self.summary = rawstream.summarize(ev)
        log(f"  {self.summary['photons']} photon records, {self.summary['sync_records']} sync records, "
            f"markers {self.summary['markers_by_bits(frame,line,pixel)']} ({time.time()-t0:.1f} s)")

        # --- which elements are photon channels -------------------------------
        n_el = r["n_elements"]
        is_pix = ev.ph_element < n_el
        other = np.unique(ev.ph_element[~is_pix], return_counts=True)
        self.excluded_other_elements = {int(a): int(b) for a, b in zip(*other)}
        if self.excluded_other_elements:
            r.warn(f"records from elements {self.excluded_other_elements} are neither SPAD pixels "
                   f"(0..{n_el-1}) nor the sync element; excluded (Abberior mentions 'dummy' and "
                   "'overflow' channels)")
        if ev.gen_type.size:
            r.warn(f"{ev.gen_type.size} 'general' messages present (types "
                   f"{sorted(set(ev.gen_type.tolist()))}); they are not interpreted")
        el = ev.ph_element[is_pix]
        code = ev.ph_code[is_pix]
        coarse = ev.ph_coarse[is_pix]

        # --- TDC calibration -----------------------------------------------------
        if r["tdc_calibration_file"]:
            cal = tdc.TdcCalibration.load(r["tdc_calibration_file"])
            filled = cal.fill_gaps(el, code)
            if filled:
                r.warn(f"loaded TDC table has no value for some codes used in this file; "
                       f"interpolated between neighbouring codes for {sum(filled.values())} photons "
                       f"(per element: {filled})")
            r.set("tdc_gap_filled_photons", filled, "loaded TDC table vs this file")
        else:
            cal = tdc.TdcCalibration.from_histograms(el, code, ev.sync.code_hist, r["sync_element"])
        self.cal = cal

        # --- pulse train: sync records are re-read from the file block-wise ----------
        pt = timing.streamed_pulse_train(ev.sync, cal.tables[r["sync_element"]],
                                         r["sync_divider"], r["sync_smoothing_window"])
        self.pt = pt
        if pt.n_missing_sync:
            r.warn(f"{pt.n_missing_sync} sync records missing in the stream (bridged; "
                   "in the test data they coincide with line markers)")

        # --- markers ----------------------------------------------------------------
        bits = ev.mk_bits
        self.line_t = ev.mk_coarse[(bits & 2) > 0].astype(np.float64)
        self.frame_t = ev.mk_coarse[(bits & 4) > 0].astype(np.float64)
        n_pix_mk = int(((bits & 1) > 0).sum())
        if n_pix_mk:
            r.warn(f"{n_pix_mk} pixel markers present; they are kept in the report but pixel "
                   "positions are computed from line start + dwell time")
        if self.line_t.size < 2:
            raise ConversionError("fewer than 2 line markers in the stream; cannot reconstruct an image")

        # --- laser frequency / time base ---------------------------------------------
        self._resolve_time_base(meta)

        # --- photons: absolute time incl. delay corrections, block-wise -----------
        res = r["microtime_resolution_ps"] * 1e-12
        # a remainder below 1/1000 bin (ppb-level frequency noise) gives no extra bin
        self.n_bins = int(np.ceil(self.Tp / res - 1e-3))
        if self.n_bins > 32768:
            raise ConversionError(f"microtime resolution {res*1e12} ps gives {self.n_bins} bins per "
                                  "laser period; the T3 record allows 32768. Use a coarser resolution.")
        delays_ps = np.asarray(r["element_delays_ps"], float) + r["common_delay_ps"]
        n = el.size
        pulse = np.empty(n, np.int64)
        dtime = np.empty(n, np.uint16)
        element = np.empty(n, np.uint8)
        t_all = np.empty(n, np.float64) if keep_validation_data else None
        ph_all = np.empty(n, np.float64) if keep_validation_data else None
        # scan markers are located in the same sweep as the photons (the pulse
        # train is recomputed block by block, so one sweep is much cheaper than two)
        off_ticks = r["line_start_offset_us"] * 1e-6 / self.tick_s
        mk_t = np.concatenate([self.line_t + off_ticks, self.frame_t])
        mk_order = np.argsort(mk_t, kind="stable")
        mk_pulse = np.empty(mk_t.size, np.int64)
        mk_done = 0
        j = 0
        for s0 in range(0, n, self.BLOCK):
            s1 = min(n, s0 + self.BLOCK)
            e = el[s0:s1]
            t = coarse[s0:s1] + cal.fine(e, code[s0:s1]) + delays_ps[e] * 1e-12 / self.tick_s
            keep = ~pt.outside_mask(t)
            t, e = t[keep], e[keep]
            # markers up to the end of this block's time range go along
            mk_new = mk_done + int(np.searchsorted(mk_t[mk_order[mk_done:]], t.max(), side="right")) \
                if t.size else mk_done
            mk_idx = mk_order[mk_done:mk_new]
            k = t.size
            p, phase = pt.locate(np.concatenate([t, mk_t[mk_idx]]))
            mk_pulse[mk_idx] = p[k:]
            p, phase = p[:k], phase[:k]
            mk_done = mk_new
            pulse[j:j + k] = p
            dtime[j:j + k] = np.minimum(np.floor(phase * self.Tp / res), self.n_bins - 1)
            element[j:j + k] = e
            if keep_validation_data:
                t_all[j:j + k] = t
                ph_all[j:j + k] = phase
            j += k
        if mk_done < mk_t.size:
            mk_idx = mk_order[mk_done:]
            mk_pulse[mk_idx], _ = pt.locate(mk_t[mk_idx])
        del el, code, coarse, ev
        self.n_dropped_outside_sync = int(n - j)
        if self.n_dropped_outside_sync:
            log(f"  {self.n_dropped_outside_sync} photons before the first / after the last sync "
                "record are dropped (no pulse reference)")
            pulse, dtime, element = pulse[:j].copy(), dtime[:j].copy(), element[:j].copy()
            if keep_validation_data:
                t_all, ph_all = t_all[:j].copy(), ph_all[:j].copy()
        # sort by laser pulse once (stable: equal pulses keep the stream order),
        # so that each output file only needs a merge with the few markers
        if pulse.size and not _is_sorted(pulse):
            order = np.argsort(pulse, kind="stable")
            pulse = pulse[order]
            dtime = dtime[order]
            element = element[order]
            if keep_validation_data:
                t_all = t_all[order]
                ph_all = ph_all[order]
            del order
        self.pulse, self.dtime, self.element = pulse, dtime, element
        self.t, self.phase = t_all, ph_all

        # --- markers in pulse units ---------------------------------------------------
        self.line_pulse = mk_pulse[:self.line_t.size]
        self.frame_pulse = mk_pulse[self.line_t.size:]
        self.line_len_pulses = int(round(r["pixels_x"] * r["pixel_dwell_us"] * 1e-6 / self.Tp))
        self.stop_pulse = self.line_pulse + self.line_len_pulses
        pt.release()            # the pulse train is no longer needed; free its memory

        # line pass index and sequence for every photon (validation / split)
        self.pass_of_photon = np.empty(self.pulse.size, np.int32)
        for s0 in range(0, self.pulse.size, self.BLOCK):
            self.pass_of_photon[s0:s0 + self.BLOCK] = np.searchsorted(
                self.line_pulse, self.pulse[s0:s0 + self.BLOCK], side="right") - 1
        self.seq_of_pass = self._sequence_pattern(self.line_pulse.size)
        self._check_scan_structure()

    # ------------------------------------------------------------------------- #
    def _resolve_time_base(self, meta):
        r, pt = self.r, self.pt
        tpp = pt.ticks_per_pulse()
        # mean line-marker interval, excluding frame gaps (> 1% off the median);
        # the mean avoids the 1-tick quantisation of individual intervals
        iv = np.diff(self.line_t)
        iv = iv[np.abs(iv / np.median(iv) - 1) < 0.01]
        pulses_per_line = float(iv.mean()) / tpp
        f_line = r["line_frequency_hz_msr"]
        f_user = r["laser_frequency_hz"]
        f_from_lines = pulses_per_line * f_line if f_line else None
        if f_user not in (AUTO, None):
            f = float(f_user)
            src = "config"
            if f_from_lines and abs(f / f_from_lines - 1) > 2e-6:
                r.warn(f"configured laser frequency {f:.1f} Hz differs from the value derived from "
                       f"the scanner line period ({f_from_lines:.1f} Hz)")
        elif f_from_lines:
            f = f_from_lines
            src = ("derived: (laser pulses per line-marker interval, from the sync channel) x "
                   "(line frequency from .msr)")
        else:
            f = 40e6
            src = "DEFAULT nominal 40 MHz (no .msr, no config value) - verify!"
            r.warn("laser frequency not known; using nominal 40 MHz. Absolute times/pixel dwell "
                   "may be off by a few ppm; set timing.laser_frequency_hz")
        r.set("laser_frequency_hz", f, src)
        if meta:
            nominal = meta.get("sync_frequency_hz")
            ext = meta.get("ext_frequencies_hz") or []
            inp = meta.get("sync_source_input")
            ms = meta.get("metadata_label", ".msr metadata")
            r.set("msr_sync_source_input", inp, ms)
            r.set("msr_ext_frequencies_hz", ext, ms)
            if inp not in (None, "-1") and ext:
                # external sync: Imspector's measured frequency of that input
                try:
                    measured = ext[int(inp)]
                    if abs(measured / f - 1) > 2e-6:
                        r.warn(f"laser frequency {f:.1f} Hz vs Imspector-measured external sync "
                               f"{measured:.1f} Hz differ by {1e6*(measured/f-1):.2f} ppm")
                except (ValueError, IndexError):
                    pass
            elif nominal and abs(nominal / f - 1) > 2e-6:
                r.warn(f"laser frequency {f:.1f} Hz differs from the internal sync setting "
                       f"{nominal:.0f} Hz by {1e6*(nominal/f-1):.2f} ppm")
        self.Tp = 1.0 / f
        self.tick_s = self.Tp / tpp
        r.set("coarse_tick_ns", self.tick_s * 1e9, "measured (sync channel + laser frequency)")
        r.set("pulses_per_line_marker_interval", pulses_per_line, "measured")
        r.set("sync_jitter_ps_rms", pt.jitter_ticks_rms * self.tick_s * 1e12,
              "measured (per sync record, before smoothing)")

    def _sequence_pattern(self, n_passes):
        acc = [s["line_accumulation"] for s in self.r["sequences"]]
        pat = np.repeat(np.arange(len(acc)), acc)
        reps = -(-n_passes // pat.size)
        return np.tile(pat, reps)[:n_passes]

    def _check_scan_structure(self):
        r = self.r
        acc_total = sum(s["line_accumulation"] for s in r["sequences"])
        expected = r["pixels_y"] * acc_total * r["frames"]
        n = self.line_t.size
        if n != expected:
            msg = (f"stream has {n} line markers, but pixels_y x sum(line_accumulation) x frames = "
                   f"{r['pixels_y']} x {acc_total} x {r['frames']} = {expected}")
            if n < expected:
                r.warn(msg + " -> acquisition incomplete; missing lines stay empty")
            else:
                raise ConversionError(msg + ". Check pixels_y, sequences/line_accumulation and frames.")
        if self.frame_t.size != r["frames"]:
            r.warn(f"stream has {self.frame_t.size} frame markers, config frames = {r['frames']}")
        intervals = np.diff(self.line_pulse)
        if intervals.size and intervals.min() < self.line_len_pulses:
            raise ConversionError(
                f"line markers are {intervals.min()*self.Tp*1e3:.3f} ms apart, shorter than "
                f"pixels_x x dwell = {self.line_len_pulses*self.Tp*1e3:.3f} ms; dwell or pixels_x is wrong")
        flyback = (np.median(intervals) - self.line_len_pulses) * self.Tp
        r.set("flyback_ms", flyback * 1e3, "measured (line period - pixels_x x dwell)")
        # photon-rate check: the rate after the line window should drop (flyback)
        n_inside = n_after = 0
        med_iv = np.median(intervals)
        for s0 in range(0, self.pulse.size, self.BLOCK):
            k = np.clip(self.pass_of_photon[s0:s0 + self.BLOCK], 0, None)
            rel = self.pulse[s0:s0 + self.BLOCK] - self.line_pulse[k]
            n_inside += int(((rel >= 0) & (rel < self.line_len_pulses)).sum())
            n_after += int(((rel >= self.line_len_pulses) & (rel < med_iv)).sum())
        if n_after and flyback > 0:
            rate_in = n_inside / (self.line_len_pulses)
            rate_after = n_after / (med_iv - self.line_len_pulses)
            r.set("flyback_rate_ratio", float(rate_after / rate_in), "measured")
            if rate_after > 0.5 * rate_in:
                r.warn(f"photon rate after the line window is {rate_after/rate_in:.0%} of the rate "
                       "inside it; line markers may not mark the line start, or dwell is wrong")


# --------------------------------------------------------------------------- #
# PTU assembly
# --------------------------------------------------------------------------- #
def _tags(acq: DecodedAcquisition, mode: str, n_records_placeholder: int, n_channels: int,
          pix_y: int, seq_desc: str, source: str, comment: str, hw_type: str):
    r = acq.r
    acq_time_ms = float((acq.pulse.max() if acq.pulse.size else 0) * acq.Tp * 1e3)
    tags = [
        ("File_GUID", new_guid()),
        ("File_CreatingTime", dt.datetime.now()),
        ("File_Comment", comment),
        ("CreatorSW_Name", "matrix2ptu"),
        ("CreatorSW_Version", __version__),
        ("Measurement_Mode", 3),
        ("Measurement_SubMode", 3),
        ("HW_Type", hw_type),
        ("HW_InpChannels", n_channels + 1),
        ("HW_Markers", 4),
        ("MeasDesc_BinningFactor", 1),
        ("MeasDesc_Offset", 0),
        ("MeasDesc_GlobalResolution", float(acq.Tp)),
        ("MeasDesc_Resolution", float(r["microtime_resolution_ps"] * 1e-12)),
        ("MeasDesc_AcquisitionTime", int(round(acq_time_ms))),
        ("TTResult_SyncRate", int(round(1.0 / acq.Tp))),
        ("TTResult_StopReason", 0),
        ("TTResult_NumberOfRecords", int(n_records_placeholder)),
        ("TTResultFormat_TTTRRecType", REC_TYPE_HYDRAHARP2_T3),
        ("TTResultFormat_BitsPerRecord", 32),
        ("TTResultFormat_WrapAround", WRAP),
        ("ImgHdr_Dimensions", 3),
        ("ImgHdr_Ident", 13),            # same scanner identifier Imspector writes
        ("ImgHdr_PixX", int(r["pixels_x"])),
        ("ImgHdr_PixY", int(pix_y)),
        ("ImgHdr_PixResol", float(r["pixel_size_um"])),
        ("ImgHdr_TimePerPixel", float(r["pixel_dwell_us"] * 1e-3)),   # ms
        ("ImgHdr_X0", float(r["x0_um"])),
        ("ImgHdr_Y0", float(r["y0_um"])),
        ("ImgHdr_LineStart", 1),
        ("ImgHdr_LineStop", 2),
        ("ImgHdr_Frame", 3),
        ("ImgHdr_BiDirect", False),
        ("ImgHdr_SinCorrection", 0),
        # --- provenance (custom tags; ignored by readers that do not know them)
        ("UsrMX_SourceFile", os.path.basename(source)),
        ("UsrMX_OutputMode", mode),
        ("UsrMX_ChannelMeaning",
         "detector element 0..22" if mode == "detector_resolved" else "all elements pooled into channel 0"),
        ("UsrMX_Sequences", seq_desc),
        ("UsrMX_LineOrder", "line_interleaved: per y-line, sequences in order, each repeated line_accumulation times"),
        ("UsrMX_DwellConvention", "ImgHdr_TimePerPixel is per line pass"),
        ("UsrMX_LineStopSynthesized", True),
        ("UsrMX_TdcCalibration", acq.cal.method),
        ("UsrMX_LaserFreqSource", str(r.source["laser_frequency_hz"])[:200]),
        ("UsrMX_ElementDelaysApplied", bool(r["apply_element_delays"])),
        ("UsrMX_CommonDelay_ps", float(r["common_delay_ps"])),
        ("UsrMX_SyncJitter_ps", float(r["sync_jitter_ps_rms"])),
    ]
    for i, d in enumerate(r["element_delays_ps"]):
        tags.append(("UsrMX_ElementDelay_ps", float(d), i))
    return tags


def _is_sorted(a: np.ndarray, block: int = 4_000_000) -> bool:
    for s0 in range(0, a.size - 1, block):
        if (np.diff(a[s0:s0 + block + 1]) < 0).any():
            return False
    return True


def _assemble(acq: DecodedAcquisition, photon_mask: np.ndarray | None, pooled: bool,
              pass_mask: np.ndarray | None):
    """Merge photons and markers into time-ordered arrays: by nsync, markers
    before photons at equal nsync, otherwise in input order (photons are
    already sorted by pulse, so only the few markers need to be inserted)."""
    sel = (lambda a: a) if photon_mask is None else (lambda a: a[photon_mask])
    ph_n = sel(acq.pulse)
    ph_d = sel(acq.dtime)
    ls = acq.line_pulse if pass_mask is None else acq.line_pulse[pass_mask]
    le = acq.stop_pulse if pass_mask is None else acq.stop_pulse[pass_mask]
    fr = acq.frame_pulse
    mk_n = np.concatenate([fr, ls, le]).astype(np.int64)
    mk_c = np.concatenate([np.full(fr.size, MARK_FRAME), np.full(ls.size, MARK_LINE_START),
                           np.full(le.size, MARK_LINE_STOP)]).astype(np.uint8)
    o = np.argsort(mk_n, kind="stable")
    mk_n, mk_c = mk_n[o], mk_c[o]
    pos = np.searchsorted(ph_n, mk_n, side="left")
    n = np.insert(ph_n, pos, mk_n)
    d = np.insert(ph_d, pos, np.zeros(mk_n.size, np.uint16))
    s = np.zeros(n.size, np.uint8)
    s[pos + np.arange(pos.size)] = 1
    if pooled:
        c = np.zeros(n.size, np.uint8)
        c[s == 1] = mk_c
    else:
        c = np.insert(sel(acq.element), pos, mk_c)
    return n, c, d, s


def convert(config_path: str, only_validate_inputs: bool = False,
            keep_validation_data: bool = False) -> dict:
    """Run the conversion.  keep_validation_data=True keeps the per-photon
    arrays that validate.validate() needs (16 bytes/photon more memory)."""
    cfg = load_config(config_path)
    raw = source_path(cfg)
    kind = "obf" if cfg["input"]["obf"] else "raw"
    if not os.path.exists(raw):
        raise ConversionError(f"{kind} file not found: {raw}")
    source_info = {"path": raw, "container": kind}
    meta = None
    meta_diffs = []
    if kind == "obf":
        # a streamed .obf carries the event words AND the Imspector parameters of
        # the acquisition; those are the primary metadata (an .msr is optional)
        from .obf import find_event_stream
        es = find_event_stream(raw)         # clear error if this is an image-only OBF
        meta = read_msr_metadata(raw)
        meta["metadata_label"] = "embedded .obf parameters"
        source_info.update({"obf_stack": es.stack_name, "obf_byte_offset": es.byte_offset,
                            "obf_stack_words": es.n_words,
                            "matrix_streaming_format": meta.get("streaming_format")})
    if cfg["input"]["msr"]:
        if not os.path.exists(cfg["input"]["msr"]):
            raise ConversionError(f"msr file not found: {cfg['input']['msr']}")
        msr_meta = read_msr_metadata(cfg["input"]["msr"])
        if meta is None:
            meta = msr_meta
        else:
            meta_diffs = compare_metadata(meta, msr_meta)
    r = resolve(cfg, meta)
    for d in meta_diffs:
        r.warn(f"embedded .obf parameters and .msr disagree: {d}  (the .msr may belong to "
               "another acquisition; the .obf values are used)")
    if kind == "obf" and cfg["input"]["msr"] and not meta_diffs:
        r.set("msr_vs_obf_parameters", "agree", "cross-check")
    acq = DecodedAcquisition(cfg, r, meta, keep_validation_data)

    out_dir = cfg["output"]["directory"]
    os.makedirs(out_dir, exist_ok=True)
    base = cfg["output"]["basename"] or os.path.splitext(os.path.basename(raw))[0]
    seqs = r["sequences"]
    seq_desc = ";".join(f"{i}:{s['name']}:x{s['line_accumulation']}" for i, s in enumerate(seqs))

    files = []
    if not only_validate_inputs:
        for mode in cfg["output"]["modes"]:
            if mode not in ("detector_resolved", "pooled"):
                raise ConversionError(f"unknown output mode '{mode}'")
            pooled = mode == "pooled"
            n_ch = 1 if pooled else int(acq.element.max()) + 1
            targets = [(None, None, None, r["pixels_y"], seq_desc, "")]
            if cfg["output"]["split_sequences"] and len(seqs) > 1:
                targets = []
                for i, s in enumerate(seqs):
                    pm = acq.seq_of_pass == i
                    phm = np.zeros(acq.pulse.size, bool)
                    for s0 in range(0, phm.size, acq.BLOCK):
                        k = acq.pass_of_photon[s0:s0 + acq.BLOCK]
                        phm[s0:s0 + acq.BLOCK] = (k >= 0) & pm[np.clip(k, 0, None)]
                    targets.append((i, phm, pm, r["pixels_y"] * s["line_accumulation"],
                                    f"0:{s['name']}:x{s['line_accumulation']}",
                                    f"_seq{i}_{_safe(s['name'])}"))
            for seq_i, phm, pm, pix_y, sdesc, suffix in targets:
                n, c, d, s = _assemble(acq, phm, pooled, pm)
                fname = os.path.join(out_dir, f"{base}_{mode}{suffix}.ptu")
                if os.path.exists(fname) and not cfg["output"]["overwrite"]:
                    raise ConversionError(f"{fname} exists; set output.overwrite: true")
                comment = (f"Converted from {os.path.basename(raw)} by matrix2ptu {__version__}. "
                           f"Mode {mode}. Sequences {sdesc}. LineStop markers synthesised at "
                           f"line start + {r['pixels_x']} x {r['pixel_dwell_us']} us.")
                tags = _tags(acq, mode, 0, n_ch, pix_y, sdesc, raw, comment,
                             cfg["output"]["hw_type"])     # record count patched on write
                if kind == "obf":
                    tags.insert(tags.index(next(t for t in tags if t[0] == "UsrMX_OutputMode")),
                                ("UsrMX_SourceFormat",
                                 f"streamed OBF, event stack '{source_info['obf_stack']}' "
                                 f"(Matrix1 streaming format {source_info['matrix_streaming_format']})"))
                n_rec = write_ptu_t3(fname, tags, n, c, d, s)
                n_ph = int(acq.pulse.size if phm is None else phm.sum())
                if pooled:
                    per_ch = {0: n_ph} if n_ph else {}
                else:
                    cnt = np.bincount(acq.element if phm is None else acq.element[phm])
                    per_ch = {int(k): int(v) for k, v in enumerate(cnt) if v}
                info = {
                    "file": fname, "mode": mode, "sequence": seq_i,
                    "records": int(n_rec), "photons": n_ph,
                    "photons_per_channel": per_ch,
                    "line_start_markers": int(acq.line_pulse.size if pm is None else pm.sum()),
                    "line_stop_markers": int(acq.line_pulse.size if pm is None else pm.sum()),
                    "frame_markers": int(acq.frame_pulse.size),
                    "overflow_records": int(n_rec - n.size),
                    "img_pix_y": int(pix_y),
                }
                files.append(info)
                del n, c, d, s
                log(f"wrote {fname}: {n_rec} records ({info['photons']} photons)")
        acq.cal.save(os.path.join(out_dir, f"{base}_tdc_calibration.npz"))

    report = {
        "matrix2ptu_version": __version__,
        "config_file": config_path,
        "source": source_info,
        "raw_summary": acq.summary,
        "parameters": {k: {"value": _jsonable(v), "source": r.source.get(k)} for k, v in r.values.items()},
        "derived": {
            "laser_period_ns": acq.Tp * 1e9,
            "microtime_bins_per_period": acq.n_bins,
            "line_length_pulses": acq.line_len_pulses,
            "n_line_markers": int(acq.line_t.size),
            "n_frame_markers": int(acq.frame_t.size),
            "photons_dropped_outside_sync_range": acq.n_dropped_outside_sync,
            "records_excluded_other_elements": acq.excluded_other_elements,
            "missing_sync_records_bridged": acq.pt.n_missing_sync,
        },
        "warnings": r.warnings,
        "outputs": files,
    }
    with open(os.path.join(out_dir, f"{base}_conversion_report.json"), "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2)
    _print_report(report)
    return {"report": report, "acq": acq, "cfg": cfg, "meta": meta}


def _safe(s):
    return "".join(ch if ch.isalnum() else "_" for ch in s)


def _jsonable(v):
    if isinstance(v, (np.floating, np.integer)):
        return v.item()
    return v


def _print_report(rep):
    log("parameters used:")
    for k, v in rep["parameters"].items():
        val = v["value"]
        if isinstance(val, list) and len(val) > 6:
            val = f"[{len(val)} values]"
        log(f"  {k:32s} = {val}   ({v['source']})")
    for k, v in rep["derived"].items():
        log(f"  {k:32s} = {v}")
    if rep["warnings"]:
        log("WARNINGS:")
        for w in rep["warnings"]:
            log(f"  ! {w}")
