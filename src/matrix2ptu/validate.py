"""Validation of converted PTU files against the source data.

Independent readers are used on purpose:
  * ptufile (C. Gohlke)  - header parsing and record decoding
  * tttrkit              - the PTU reader + image reconstructor used by FLOPA
Neither shares code with matrix2ptu's writer.

The PTU files are decoded block by block (two passes per file), so memory
does not grow with the file size beyond the converter's own per-photon arrays.

Checks (all results go into a JSON report and are printed):
  S1  structure: header tags readable, record type, record count, nsync
      monotonic after overflow correction, dtime < bins/period
  S2  event counts: photons per channel, markers by type vs. the source
  S3  event round trip: (nsync, dtime, channel) of every photon equals what the
      converter computed from the raw stream (exact, and in the same order)
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


CHUNK = 2_000_000        # PTU records / photons per block (bounds temporary memory)
T3_WRAP = 1024           # HydraHarp-v2 / MultiHarp T3: 10-bit nsync field


def iter_ptu_records(path, chunk: int = CHUNK):
    """Decode a PTU with ptufile, block by block.

    Yields (index of the first record, records) with absolute `time`.  ptufile
    decodes each block starting from time 0, so the overflow count of the
    previous blocks is added.  Each block is cut after its last photon/marker
    record (the following overflow records start the next block); the carry is
    then that record's absolute time rounded down to a multiple of the T3 wrap,
    exact because a photon/marker time is wraps x 1024 + its 10-bit nsync.
    """
    import ptufile
    p = ptufile.PtuFile(path)
    try:
        raw = p.read_records(copy=False)
        n = raw.size
        s0, base = 0, 0
        while s0 < n:
            s1 = min(n, s0 + chunk)
            while True:          # a block must contain a photon/marker (or end the file)
                d = p.decode_records(raw[s0:s1])
                ev = np.flatnonzero((d["channel"] >= 0) | (d["marker"] > 0))
                if ev.size or s1 == n:
                    break
                s1 = min(n, s1 + chunk)
            if s1 < n and ev[-1] + 1 < d.size:
                d = d[:ev[-1] + 1]                  # trailing overflows go to the next block
            t = d["time"].astype(np.int64) + base
            if ev.size:
                base = int(t[ev[-1]]) // T3_WRAP * T3_WRAP
            yield s0, d, t
            s0 += d.size
    finally:
        p.close()


def ptu_header(path):
    import ptufile
    p = ptufile.PtuFile(path)
    try:
        return dict(p.tags), int(p.record_type)
    finally:
        p.close()


def _image_add(img, tp, ls, le, pixels, seq_of_pass, passes_per_line, n_lines):
    """Add photons at nsync times tp to img[seq, y, x] using the PTU's markers:
    LineStart/LineStop delimit each line pass; pixel = floor(phase * pixels)
    with phase measured in nsync (as FLOPA/SymPhoTime do).  Line pass k
    belongs to y = k // passes_per_line and sequence seq_of_pass[k]."""
    k = np.searchsorted(ls, tp, side="right") - 1
    kk = np.clip(k, 0, None)
    ok = (k >= 0) & (tp < le[kk])
    x = np.floor((tp - ls[kk]) / (le[kk] - ls[kk]) * pixels).astype(np.int64)
    ok &= (x >= 0) & (x < pixels)
    y = kk // passes_per_line
    s = seq_of_pass[np.minimum(kk, seq_of_pass.size - 1)]
    ok &= y < n_lines
    np.add.at(img, (s[ok], y[ok], x[ok]), 1)


def image_from_raw(acq: DecodedAcquisition, use_pulses=False):
    """Image directly from raw time stamps (no PTU): line marker + k*dwell.
    Computed in blocks of photons."""
    r = acq.r
    nseq = len(r["sequences"])
    unit = sum(s["line_accumulation"] for s in r["sequences"])
    line_t = acq.line_t + r["line_start_offset_us"] * 1e-6 / acq.tick_s
    d_ps = np.asarray(r["element_delays_ps"]) + r["common_delay_ps"]
    img = np.zeros((nseq, r["pixels_y"], r["pixels_x"]), np.int64)
    for s0 in range(0, acq.pulse.size, CHUNK):
        sl = slice(s0, s0 + CHUNK)
        if use_pulses:
            k = acq.pass_of_photon[sl]
            rel = (acq.pulse[sl] - acq.line_pulse[np.clip(k, 0, None)]) / acq.line_len_pulses
        else:
            # use the photon's raw time stamp *without* delay correction
            t_raw = acq.t[sl] - d_ps[acq.element[sl]] * 1e-12 / acq.tick_s
            k = np.searchsorted(line_t, t_raw, side="right") - 1
            rel = (t_raw - line_t[np.clip(k, 0, None)]) * acq.tick_s / (r["pixels_x"] * r["pixel_dwell_us"] * 1e-6)
        x = np.floor(rel * r["pixels_x"]).astype(np.int64)
        ok = (k >= 0) & (x >= 0) & (x < r["pixels_x"])
        kk = np.clip(k, 0, None)
        y = kk // unit
        s = acq.seq_of_pass[np.minimum(kk, acq.seq_of_pass.size - 1)]
        ok &= y < r["pixels_y"]
        np.add.at(img, (s[ok], y[ok], x[ok]), 1)
    return img


def _roundtrip_sorted(path, acq, exp_ch):
    """Order-independent round trip (the pre-streaming S3; needs the whole file
    in memory).  Only used when the in-order comparison fails, to tell an
    ordering difference from a content difference."""
    import ptufile
    p = ptufile.PtuFile(path)
    rec = p.decode_records()
    p.close()
    t = rec["time"].astype(np.int64)
    ch = rec["channel"]
    is_ph = (ch >= 0) & (rec["marker"] == 0)
    a = np.lexsort((acq.dtime, exp_ch, acq.pulse))
    b = np.lexsort((rec["dtime"][is_ph], ch[is_ph], t[is_ph]))
    return (np.array_equal(acq.pulse[a], t[is_ph][b]) and
            np.array_equal(acq.dtime[a].astype(np.int64), rec["dtime"][is_ph][b].astype(np.int64)) and
            np.array_equal(exp_ch[a].astype(np.int64), ch[is_ph][b].astype(np.int64)))


class _Moments:
    """Count, mean and sum of squared deviations, combined block by block
    (Chan et al.), for the std of values that are never all in memory."""

    def __init__(self):
        self.n, self.mean, self.m2 = 0, 0.0, 0.0

    def add(self, x):
        if x.size == 0:
            return
        n2, mu2 = x.size, float(x.mean())
        m2 = float(((x - mu2) ** 2).sum())
        n = self.n + n2
        delta = mu2 - self.mean
        self.mean += delta * n2 / n
        self.m2 += m2 + delta ** 2 * self.n * n2 / n
        self.n = n

    def std(self):
        return float(np.sqrt(self.m2 / self.n)) if self.n else float("nan")


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
    if acq.t is None:
        raise ValueError("validation needs the per-photon data kept during conversion: "
                         "call convert(..., keep_validation_data=True)")
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
        tags, rtype = ptu_header(path)
        combined = f["sequence"] is None
        if f["mode"] == "detector_resolved":
            exp_ch = acq.element
        else:
            exp_ch = None                     # pooled: every photon in channel 0
        res = tags["MeasDesc_Resolution"]
        Tg = tags["MeasDesc_GlobalResolution"]

        # ---- pass 1: structure, counts, round trip (S1-S4); collect markers
        n_rec = 0
        mono = True
        last_ev = None
        dmax = -1
        ch_counts = np.zeros(0, np.int64)
        m_counts = {1: 0, 2: 0, 4: 0}
        n_ovf = 0
        n_ph = 0
        mk_t = {1: [], 2: [], 4: []}
        same = True
        dres_max = 0.0
        dres_m = _Moments()
        for _, d, t in iter_ptu_records(path):
            n_rec += d.size
            ch, mk = d["channel"], d["marker"]
            # overflow records carry no event time (ptufile reports base+count for
            # them), so monotonicity is checked on photons and markers only
            te = t[(ch >= 0) | (mk > 0)]
            if te.size:
                mono &= bool((np.diff(te) >= 0).all()) and (last_ev is None or te[0] >= last_ev)
                last_ev = te[-1]
            is_ph = (ch >= 0) & (mk == 0)
            k = int(is_ph.sum())
            tp = t[is_ph]
            dt = d["dtime"][is_ph]
            chp = ch[is_ph]
            if k:
                dmax = max(dmax, int(dt.max()))
            bc = np.bincount(chp.astype(np.int64))
            if bc.size > ch_counts.size:
                bc[:ch_counts.size] += ch_counts
                ch_counts = bc
            else:
                ch_counts[:bc.size] += bc
            for bit in (1, 2, 4):
                m = (mk & bit) > 0
                m_counts[bit] += int(m.sum())
                mk_t[bit].append(t[m])
            n_ovf += int(((ch < 0) & (mk == 0)).sum())
            if combined and same:
                # S3: photons are written in the converter's order (sorted by
                # pulse, ties in stream order), so compare in file order
                sl = slice(n_ph, n_ph + k)
                same = (n_ph + k <= acq.pulse.size and
                        np.array_equal(acq.pulse[sl], tp) and
                        np.array_equal(acq.dtime[sl].astype(np.int64), dt.astype(np.int64)) and
                        (np.array_equal(exp_ch[sl].astype(np.int64), chp.astype(np.int64))
                         if exp_ch is not None else bool((chp == 0).all())))
                if same:
                    # S4: nsync equality is proven exactly in S3; compare the micro
                    # time part (avoids float rounding of absolute times of ~100 s)
                    dres = ((dt + 0.5) * res - acq.phase[sl] * acq.Tp) * 1e12
                    if dres.size:
                        dres_max = max(dres_max, float(np.abs(dres).max()))
                    dres_m.add(dres)
            n_ph += k
        mk_t = {b: np.concatenate(v) if v else np.empty(0, np.int64) for b, v in mk_t.items()}

        check("S1 structure (ptufile)", rtype == 0x01010304 and mono and
              tags["TTResult_NumberOfRecords"] == f["records"] and n_rec == f["records"],
              f"record type {rtype:#x}, {tags['TTResult_NumberOfRecords']} records "
              f"(written {f['records']}), time monotonic={mono}, "
              f"SyncRate {tags['TTResult_SyncRate']} Hz, resolution {tags['MeasDesc_Resolution']*1e12:g} ps")
        check("S1 dtime range", dmax < acq.n_bins, f"max dtime {dmax} < bins per period {acq.n_bins}")
        # ---- S2 counts
        pc = {int(i): int(v) for i, v in enumerate(ch_counts) if v}
        check("S2 photons per channel", pc == f["photons_per_channel"],
              f"PTU {sum(pc.values())} photons in {len(pc)} channels; source/converter "
              f"{f['photons']}; raw file has {acq.summary['photons']} photon records, "
              f"{rep['derived']['photons_dropped_outside_sync_range']} dropped (outside sync range)")
        check("S2 markers", m_counts == {1: f["line_start_markers"], 2: f["line_stop_markers"],
                                         4: f["frame_markers"]},
              f"PTU LineStart {m_counts[1]}, LineStop {m_counts[2]}, Frame {m_counts[4]}; "
              f"raw line markers {rep['derived']['n_line_markers']}, frame {rep['derived']['n_frame_markers']}; "
              f"non-event (overflow) records {tags['TTResult_NumberOfRecords'] - n_ph - sum(m_counts.values())}")

        if not combined:
            out["notes"].append(f"{path}: per-sequence file, image checks done on combined files only")
            continue
        # ---- S3 exact round trip
        same = same and n_ph == acq.pulse.size
        if same:
            detail = ("every photon's (nsync, dtime, channel) decoded by ptufile equals the "
                      "converter's values")
        else:
            full_exp = exp_ch if exp_ch is not None else np.zeros_like(acq.element)
            if n_ph == acq.pulse.size and _roundtrip_sorted(path, acq, full_exp):
                detail = "same photons, but in a different order than the converter's (sorted by pulse)"
            else:
                detail = "MISMATCH"
        check("S3 event round trip", same, detail)
        # ---- S4 timing preservation
        if same:
            check("S4 timing preservation", dres_max <= res * 1e12 * 0.501,
                  f"PTU micro time (bin centre) - raw-derived micro time: max |d| {dres_max:.2f} ps, "
                  f"rms {dres_m.std():.2f} ps (bin {res*1e12:g} ps; expected limit = half a bin)")
        else:
            check("S4 timing preservation", False, "not evaluated: S3 failed")

        # ---- pass 2: images and histograms (I1, D1, D2)
        seq_pattern = acq.seq_of_pass
        ls, le = mk_t[1], mk_t[2]
        if ls.size != le.size:
            raise ValueError(f"{ls.size} LineStart vs {le.size} LineStop markers")
        names = list(cubes)
        res_ps = res * 1e12
        img_ptu = np.zeros((nseq, r["pixels_y"], r["pixels_x"]), np.int64)
        dec = {}                                       # D1: sequence -> 100 ps histogram
        nb = int(np.ceil(Tg / res - 1e-6))
        L = max(nb, dmax + 1)
        n_chan = ch_counts.size
        h2 = np.zeros(n_chan * L, np.int64)            # D2: channel x dtime
        for _, d, t in iter_ptu_records(path):
            is_ph = (d["channel"] >= 0) & (d["marker"] == 0)
            tp = t[is_ph]
            _image_add(img_ptu, tp, ls, le, r["pixels_x"], seq_pattern, unit, r["pixels_y"])
            dt_ph = d["dtime"][is_ph].astype(np.int64)
            if names:
                k = np.searchsorted(ls, tp, side="right") - 1
                kk = np.clip(k, 0, None)
                sq = seq_pattern[np.minimum(kk, seq_pattern.size - 1)]
                for i in range(min(nseq, len(names))):
                    sel = (k >= 0) & (sq == i)
                    hb = np.bincount((dt_ph[sel] * res_ps // 100).astype(np.int64))
                    prev = dec.get(i, np.zeros(0, np.int64))
                    if hb.size > prev.size:
                        hb[:prev.size] += prev
                        dec[i] = hb
                    else:
                        prev[:hb.size] += hb
                        dec[i] = prev
            if f["mode"] == "detector_resolved" and n_chan:
                h2 += np.bincount(d["channel"][is_ph].astype(np.int64) * L + dt_ph, minlength=h2.size)

        # ---- I1 own reconstruction from PTU vs raw
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
                del ds, pcf, flopa_img
            except Exception as e:  # report, do not hide
                check("I2 FLOPA/tttrkit reconstruction", False, f"exception: {e!r}")
        # ---- I3 vs Imspector cubes
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
            # ---- D1 decay (photons of this sequence, after the first line start)
            fb = int(round(ref.shape[-1]))
            h = np.zeros(fb, np.int64)
            hb = dec.get(i, np.zeros(0, np.int64))[:fb]
            h[:hb.size] = hb
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
            h2 = h2.reshape(n_chan, L).astype(float)
            ref_h = h2.sum(0)
            shifts, weights = [], []
            for e in np.flatnonzero(ch_counts):
                hh = h2[e]
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
