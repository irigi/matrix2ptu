"""Validation of converted PTU files against the source data.

Independent readers are used on purpose:
  * ptufile (C. Gohlke)  - header parsing and record decoding
  * tttrkit              - the PTU reader + image reconstructor used by FLOPA
Neither shares code with matrix2ptu's writer.

Checks (all results go into a JSON report and are printed):
  S1  structure: header tags readable, record type, record count, nsync
      monotonic after overflow correction, dtime < bins/period
  S2  event counts: photons per channel, markers by type vs. the source
  S3  event round trip: (nsync, dtime, channel) of every photon equals what the
      converter computed from the raw stream (exact)
  S4  timing preservation: PTU time vs. raw-derived time (residual <= 1 bin)
  I1  image from PTU (own simple reconstruction from ptufile records + PTU markers)
      vs. image computed directly from raw time stamps (no PTU involved)
  I2  image from PTU via tttrkit/FLOPA ImageReconstructor vs. I1
  I3  image vs. Imspector cube(s) from the .msr (pixel-wise)
  I4  image vs. Intensity.tif (correlation; the TIFF is a colour rendering)
  D1  decay histogram (100 ps bins) vs. Imspector cube histogram
  D2  per-element alignment after delay correction
"""
from __future__ import annotations

import json
import os

import numpy as np

from .converter import DecodedAcquisition, log


def _corr(a, b):
    a = np.asarray(a, float).ravel()
    b = np.asarray(b, float).ravel()
    return float(np.corrcoef(a, b)[0, 1])


def _read_ptufile(path):
    import ptufile
    p = ptufile.PtuFile(path)
    rec = p.decode_records()
    tags = dict(p.tags)
    rtype = int(p.record_type)
    p.close()
    return tags, rtype, rec


def image_from_ptu_records(rec, n_lines, pixels, seq_of_pass, n_seq, passes_per_line):
    """Own reconstruction from decoded PTU records using the PTU's markers:
    LineStart/LineStop delimit each line pass; pixel = floor(phase * pixels)
    with phase measured in nsync (as FLOPA/SymPhoTime do).  Line pass k
    belongs to y = k // passes_per_line and sequence seq_of_pass[k].
    Returns counts[seq, y, x]."""
    t = rec["time"].astype(np.int64)
    mk = rec["marker"]
    ch = rec["channel"]
    is_ph = (ch >= 0) & (mk == 0)
    ls = t[(mk & 1) > 0]
    le = t[(mk & 2) > 0]
    if ls.size != le.size:
        raise ValueError(f"{ls.size} LineStart vs {le.size} LineStop markers")
    tp = t[is_ph]
    k = np.searchsorted(ls, tp, side="right") - 1
    kk = np.clip(k, 0, None)
    ok = (k >= 0) & (tp < le[kk])
    x = np.floor((tp - ls[kk]) / (le[kk] - ls[kk]) * pixels).astype(np.int64)
    ok &= (x >= 0) & (x < pixels)
    y = kk // passes_per_line
    s = seq_of_pass[np.minimum(kk, seq_of_pass.size - 1)]
    ok &= y < n_lines
    img = np.zeros((n_seq, n_lines, pixels), np.int64)
    np.add.at(img, (s[ok], y[ok], x[ok]), 1)
    return img, is_ph, ok


def image_from_raw(acq: DecodedAcquisition, use_pulses=False):
    """Image directly from raw time stamps (no PTU): line marker + k*dwell."""
    r = acq.r
    nseq = len(r["sequences"])
    unit = sum(s["line_accumulation"] for s in r["sequences"])
    line_t = acq.line_t + r["line_start_offset_us"] * 1e-6 / acq.tick_s
    if use_pulses:
        k = acq.pass_of_photon
        rel = (acq.pulse - acq.line_pulse[np.clip(k, 0, None)]) / acq.line_len_pulses
    else:
        # use the photon's raw time stamp *without* delay correction
        d_ps = np.asarray(r["element_delays_ps"]) + r["common_delay_ps"]
        t_raw = acq.t - d_ps[acq.element] * 1e-12 / acq.tick_s
        k = np.searchsorted(line_t, t_raw, side="right") - 1
        rel = (t_raw - line_t[np.clip(k, 0, None)]) * acq.tick_s / (r["pixels_x"] * r["pixel_dwell_us"] * 1e-6)
    x = np.floor(rel * r["pixels_x"]).astype(np.int64)
    ok = (k >= 0) & (x >= 0) & (x < r["pixels_x"])
    kk = np.clip(k, 0, None)
    y = kk // unit
    s = acq.seq_of_pass[np.minimum(kk, acq.seq_of_pass.size - 1)]
    ok &= y < r["pixels_y"]
    img = np.zeros((nseq, r["pixels_y"], r["pixels_x"]), np.int64)
    np.add.at(img, (s[ok], y[ok], x[ok]), 1)
    return img


def flopa_reconstruct(path, r, max_detector=64):
    from tttrkit.ptuio.decoder import T3OverflowCorrector
    from tttrkit.ptuio.reader import TTTRReader
    from tttrkit.ptuio.reconstructor import ImageReconstructor, ScanConfig
    from tttrkit.ptuio.utils import estimate_tcspc_bins

    reader = TTTRReader(path)
    tags = reader.header.tags
    cfg = ScanConfig(lines=r["pixels_y"], pixels=r["pixels_x"], frames=r["frames"],
                     max_detector=max_detector,
                     line_accumulations=tuple(s["line_accumulation"] for s in r["sequences"]),
                     bidirectional=False)
    bins = estimate_tcspc_bins(tags, buffer=10)
    recon = ImageReconstructor(cfg, outputs=["photon_count", "tcspc_histogram"],
                               laser_sync_rate=tags["TTResult_SyncRate"],
                               tcspc_channels=bins,
                               tcspc_resolution=tags["MeasDesc_Resolution"])
    corr = T3OverflowCorrector(wraparound=tags.get("TTResultFormat_WrapAround", 1024))
    for chunk in reader.iter_chunks(1_000_000):
        recon.update(corr.correct(chunk))
        if recon._finished:
            break
    ds = recon.finalize()
    return ds, {"tcspc_bins": bins, "line_duration_nsync": cfg.line_duration,
                "sync_rate": tags["TTResult_SyncRate"]}


def validate(result: dict, tif_path: str | None = None, run_flopa: bool = True) -> dict:
    acq: DecodedAcquisition = result["acq"]
    rep = result["report"]
    meta = result["meta"]
    r = acq.r
    out = {"checks": [], "notes": []}

    def check(name, passed, detail):
        out["checks"].append({"check": name, "passed": bool(passed) if passed is not None else None,
                              "detail": detail})
        flag = "PASS" if passed else ("INFO" if passed is None else "FAIL")
        log(f"[{flag}] {name}: {detail}")

    nseq = len(r["sequences"])
    unit = sum(s["line_accumulation"] for s in r["sequences"])
    raw_img = image_from_raw(acq)
    raw_img_p = image_from_raw(acq, use_pulses=True)

    cubes = {}
    if meta is not None and result["cfg"]["input"]["msr"]:
        from .msr_meta import read_flim_stacks
        cubes = read_flim_stacks(result["cfg"]["input"]["msr"])

    for f in rep["outputs"]:
        path = f["file"]
        log(f"validating {path}")
        tags, rtype, rec = _read_ptufile(path)
        # ---- S1 structure
        t = rec["time"].astype(np.int64)
        # overflow records carry no event time (ptufile reports base+count for
        # them), so monotonicity is checked on photons and markers only
        is_event = (rec["channel"] >= 0) | (rec["marker"] > 0)
        mono = bool((np.diff(t[is_event]) >= 0).all())
        check("S1 structure (ptufile)", rtype == 0x01010304 and mono and
              tags["TTResult_NumberOfRecords"] == f["records"],
              f"record type {rtype:#x}, {tags['TTResult_NumberOfRecords']} records "
              f"(written {f['records']}), time monotonic={mono}, "
              f"SyncRate {tags['TTResult_SyncRate']} Hz, resolution {tags['MeasDesc_Resolution']*1e12:g} ps")
        ch = rec["channel"]
        mk = rec["marker"]
        is_ph = (ch >= 0) & (mk == 0)
        dmax = int(rec["dtime"][is_ph].max())
        check("S1 dtime range", dmax < acq.n_bins, f"max dtime {dmax} < bins per period {acq.n_bins}")
        # ---- S2 counts
        pc = {int(a): int(b) for a, b in zip(*np.unique(ch[is_ph], return_counts=True))}
        m_counts = {int(b): int((mk & b > 0).sum()) for b in (1, 2, 4)}
        n_ovf = int(((ch < 0) & (mk == 0)).sum())
        check("S2 photons per channel", pc == f["photons_per_channel"],
              f"PTU {sum(pc.values())} photons in {len(pc)} channels; source/converter "
              f"{f['photons']}; raw file has {acq.summary['photons']} photon records, "
              f"{rep['derived']['photons_dropped_outside_sync_range']} dropped (outside sync range)")
        check("S2 markers", m_counts == {1: f["line_start_markers"], 2: f["line_stop_markers"],
                                         4: f["frame_markers"]},
              f"PTU LineStart {m_counts[1]}, LineStop {m_counts[2]}, Frame {m_counts[4]}; "
              f"raw line markers {rep['derived']['n_line_markers']}, frame {rep['derived']['n_frame_markers']}; "
              f"non-event (overflow) records {tags['TTResult_NumberOfRecords'] - is_ph.sum() - sum(m_counts.values())}")

        if f["sequence"] is not None:
            out["notes"].append(f"{path}: per-sequence file, image checks done on combined files only")
            continue
        # ---- S3 exact round trip
        if f["mode"] == "detector_resolved":
            exp_ch = acq.element
        else:
            exp_ch = np.zeros_like(acq.element)
        a = np.lexsort((acq.dtime, exp_ch, acq.pulse))
        b = np.lexsort((rec["dtime"][is_ph], ch[is_ph], t[is_ph]))
        same = (np.array_equal(acq.pulse[a], t[is_ph][b]) and
                np.array_equal(acq.dtime[a].astype(np.int64), rec["dtime"][is_ph][b].astype(np.int64)) and
                np.array_equal(exp_ch[a].astype(np.int64), ch[is_ph][b].astype(np.int64)))
        check("S3 event round trip", same, "every photon's (nsync, dtime, channel) decoded by ptufile "
              "equals the converter's values" if same else "MISMATCH")
        # ---- S4 timing preservation
        res = tags["MeasDesc_Resolution"]
        Tg = tags["MeasDesc_GlobalResolution"]
        # nsync equality is proven exactly in S3; compare the micro time part
        # (avoids float rounding of absolute times of ~100 s)
        dres = ((rec["dtime"][is_ph][b] + 0.5) * res - acq.phase[a] * acq.Tp) * 1e12
        check("S4 timing preservation", np.abs(dres).max() <= res * 1e12 * 0.501,
              f"PTU micro time (bin centre) - raw-derived micro time: max |d| {np.abs(dres).max():.2f} ps, rms {dres.std():.2f} ps "
              f"(bin {res*1e12:g} ps; expected limit = half a bin)")

        # ---- I1 own reconstruction from PTU vs raw
        seq_pattern = acq.seq_of_pass
        img_ptu, _, _ = image_from_ptu_records(rec, r["pixels_y"], r["pixels_x"], seq_pattern, nseq, unit)
        same_p = np.array_equal(img_ptu, raw_img_p)
        d = np.abs(img_ptu - raw_img).sum()
        check("I1 image: PTU markers vs raw", same_p,
              f"identical to raw pulse-based assignment: {same_p}; vs raw time-stamp assignment: "
              f"sum|diff| = {d} ({d/raw_img.sum()*100:.3f}% of {raw_img.sum()} counts; boundary "
              "photons within one laser period / delay shift)")
        # ---- I2 FLOPA engine
        if run_flopa:
            try:
                ds, finfo = flopa_reconstruct(path, r)
                pcf = ds["photon_count"].values  # frame, sequence, line, pixel, channel
                flopa_img = pcf[0].sum(-1)
                same_f = np.array_equal(flopa_img, img_ptu)
                dd = np.abs(flopa_img - img_ptu).sum()
                check("I2 FLOPA/tttrkit reconstruction", dd <= 1e-3 * img_ptu.sum(),
                      f"tttrkit ImageReconstructor: {int(flopa_img.sum())} photons in image, shape "
                      f"{pcf.shape} (frame, sequence, line, pixel, channel); vs own PTU image "
                      f"identical={same_f}, sum|diff|={dd}; line_duration {finfo['line_duration_nsync']} "
                      f"nsync, tcspc bins {finfo['tcspc_bins']}")
                if f["mode"] == "detector_resolved":
                    nch = pcf.shape[-1]
                    per_ch = pcf[0].sum(axis=(0, 1, 2))
                    exp = np.array([f["photons_per_channel"].get(i, 0) for i in range(nch)])
                    out["notes"].append(f"FLOPA channels in image: {nch}; per-channel totals "
                                        f"(in-image) {per_ch.tolist()} vs all photons {exp.tolist()}")
            except Exception as e:  # report, do not hide
                check("I2 FLOPA/tttrkit reconstruction", False, f"exception: {e!r}")
        # ---- I3 vs Imspector cubes
        names = list(cubes)
        for i, s in enumerate(r["sequences"]):
            if i >= len(names):
                break
            ref = cubes[names[i]].astype(np.int64)
            R = ref.sum(-1)
            if R.shape != img_ptu[i].shape:
                check(f"I3 vs .msr '{names[i]}'", False, f"shape {R.shape} vs {img_ptu[i].shape}")
                continue
            dsum = int(np.abs(img_ptu[i] - R).sum())
            c = _corr(img_ptu[i], R)
            ident = float((img_ptu[i] == R).mean())
            passed = c > 0.99
            hint = ""
            if c < 0.5:
                hint = (" -> images are unrelated: the .msr stack is most likely NOT from the same "
                        "acquisition as the .raw file (decay comparison below is then only indicative)")
            check(f"I3 image vs .msr '{names[i]}'", passed,
                  f"PTU {int(img_ptu[i].sum())} vs Imspector {int(R.sum())} counts; identical pixels "
                  f"{ident*100:.2f}%; sum|diff| {dsum} ({dsum/max(R.sum(),1)*100:.3f}%); r = {c:.4f}" + hint)
            # ---- D1 decay
            res_ps = tags["MeasDesc_Resolution"] * 1e12
            fb = int(round(ref.shape[-1]))
            bw = 100.0 / res_ps
            sel = np.zeros(is_ph.sum(), bool)
            # photons of this sequence in the image
            dt_ph = rec["dtime"][is_ph].astype(np.int64)
            k = np.searchsorted(t[(mk & 1) > 0], t[is_ph], side="right") - 1
            kk = np.clip(k, 0, None)
            sel = (k >= 0) & (seq_pattern[np.minimum(kk, seq_pattern.size - 1)] == i)
            h = np.bincount((dt_ph[sel] * res_ps // 100).astype(np.int64), minlength=fb)[:fb]
            H = ref.sum((0, 1)).astype(float)
            hn, Hn = h / h.sum(), H / H.sum()
            l1 = float(np.abs(hn - Hn).sum())
            chi = float((((h - H * h.sum() / H.sum()) ** 2) / np.maximum(H * h.sum() / H.sum(), 1)).mean())
            check(f"D1 decay vs .msr '{names[i]}'", (l1 < 0.03) if c >= 0.5 else None,
                  f"peak bin PTU {int(h.argmax())} / Imspector {int(H.argmax())} (100 ps bins); "
                  f"normalised L1 distance {l1:.4f}; chi2/bin {chi:.2f}"
                  + ("" if c >= 0.5 else " [reference is a different acquisition: shape comparison only]"))
        # ---- I4 Intensity.tif
        if tif_path and os.path.exists(tif_path) and nseq == 1:
            import tifffile
            tif = tifffile.imread(tif_path).astype(float)
            if tif.ndim == 3:
                cs = [_corr(img_ptu[0], tif[..., j]) for j in range(tif.shape[-1])]
                check("I4 vs Intensity.tif", max(cs) > 0.9,
                      f"TIFF is RGB {tif.shape} (colour-mapped rendering, not counts); correlation of "
                      f"PTU image with R,G,B = {[round(c,3) for c in cs]}")
            else:
                check("I4 vs Intensity.tif", _corr(img_ptu[0], tif) > 0.9,
                      f"r = {_corr(img_ptu[0], tif):.4f}")
        # ---- D2 element alignment
        if f["mode"] == "detector_resolved":
            res_ps = tags["MeasDesc_Resolution"] * 1e12
            dt_ph = rec["dtime"][is_ph].astype(np.int64)
            chp = ch[is_ph]
            nb = int(np.ceil(Tg / tags["MeasDesc_Resolution"] - 1e-6))
            ref_h = np.bincount(dt_ph, minlength=nb).astype(float)
            shifts, weights = [], []
            for e in np.unique(chp):
                hh = np.bincount(dt_ph[chp == e], minlength=nb).astype(float)
                cc = np.fft.irfft(np.fft.rfft(ref_h) * np.conj(np.fft.rfft(hh)), n=nb)
                sft = int(np.argmax(cc))
                sft = sft if sft < nb / 2 else sft - nb
                shifts.append(sft * res_ps)
                weights.append(hh.sum())
            shifts = np.array(shifts)
            w = np.array(weights)
            spread = float(np.sqrt(np.average((shifts - np.average(shifts, weights=w)) ** 2, weights=w)))
            check("D2 element alignment", spread < 100,
                  f"residual per-element offsets after delay correction: weighted rms {spread:.0f} ps, "
                  f"range {shifts.min():.0f}..{shifts.max():.0f} ps (weak elements are noisier)")
            out["element_residual_offsets_ps"] = shifts.tolist()

    out_dir = result["cfg"]["output"]["directory"]
    base = os.path.splitext(os.path.basename(rep["outputs"][0]["file"]))[0].rsplit("_", 2)[0] \
        if rep["outputs"] else "validation"
    np.save(os.path.join(out_dir, "raw_image.npy"), raw_img)
    with open(os.path.join(out_dir, "validation_report.json"), "w") as fh:
        json.dump(out, fh, indent=2)
    n_fail = sum(1 for c in out["checks"] if c["passed"] is False)
    log(f"validation finished: {len(out['checks'])} checks, {n_fail} failed")
    return out
