"""Unit tests on synthetic data (no input files needed).

Run:  .venv/bin/pytest -q
"""
import os

import numpy as np
import pytest

from matrix2ptu import rawstream, tdc, timing
from matrix2ptu.ptu_writer import encode_t3, encode_t3_fast, write_ptu

ROOT = os.path.dirname(os.path.dirname(__file__))


# --------------------------------------------------------------------------- #
def _word_photon(element, code, coarse):
    return (1 << 31) | (element << 26) | (code << 16) | (coarse & 0xFFFF)


def _word_marker(frame, line, pixel, coarse):
    return (1 << 29) | (frame << 28) | (line << 27) | (pixel << 26) | (coarse & 0xFFFF)


def _synthetic_stream(tmp_path, n_sync=3000, period=78.33):
    """Sync every `period` ticks (spanning several 16-bit wraps), a few
    photons and markers, then zero padding like the real files."""
    words, truth = [], []
    t_sync = np.round(np.arange(n_sync) * period).astype(np.int64) + 100
    events = [(t, 0, ("sync",)) for t in t_sync]
    events.append((t_sync[10] + 5, 1, ("frame",)))
    events.append((t_sync[20] + 3, 1, ("line",)))
    rng = np.random.default_rng(0)
    pt = np.sort(rng.integers(t_sync[0] + 1, t_sync[-1], 500))
    for t in pt:
        events.append((int(t), 2, ("ph", int(rng.integers(0, 23)), int(rng.integers(0, 1024)))))
    events.sort(key=lambda e: (e[0], e[1]))
    for t, _, ev in events:
        if ev[0] == "sync":
            words.append(_word_photon(23, 5, t))
        elif ev[0] == "frame":
            words.append(_word_marker(1, 0, 0, t))
        elif ev[0] == "line":
            words.append(_word_marker(0, 1, 0, t))
        else:
            words.append(_word_photon(ev[1], ev[2], t))
            truth.append((t, ev[1], ev[2]))
    arr = np.array(words + [0] * 7, dtype="<u4")
    p = tmp_path / "syn.raw"
    arr.tofile(p)
    return str(p), t_sync, truth


def test_decode_unwrap(tmp_path):
    path, t_sync, truth = _synthetic_stream(tmp_path)
    ev = rawstream.decode(path)
    assert ev.n_padding_words == 7
    assert ev.sy_coarse.size == t_sync.size
    assert np.array_equal(ev.sy_coarse - ev.sy_coarse[0], t_sync - t_sync[0])
    assert t_sync[-1] > 65536 * 3            # several wraps happened
    tt = np.array([t for t, _, _ in truth])
    assert np.array_equal(ev.ph_coarse - ev.sy_coarse[0], tt - t_sync[0])
    assert np.array_equal(ev.ph_element, [e for _, e, _ in truth])
    assert np.array_equal(ev.ph_code, [c for _, _, c in truth])
    s = rawstream.summarize(ev)
    assert s["markers_by_bits(frame,line,pixel)"] == {"010": 1, "100": 1}


def test_unwrap_rejects_unordered():
    w = np.array([_word_photon(23, 0, 1000), _word_photon(23, 0, 0)], dtype=np.uint32)
    with pytest.raises(ValueError):
        rawstream.unwrap_coarse(w, max_back=10)


# --------------------------------------------------------------------------- #
def test_tdc_code_density_reversed():
    """Generate codes from a known non-uniform tap structure (reversed taps,
    4 quarters) and check that the table recovers the true positions."""
    rng = np.random.default_rng(1)
    widths = rng.uniform(0.2, 2.0, 140)
    edges = np.concatenate([[0], np.cumsum(widths)]) / widths.sum()   # within a quarter
    u = rng.uniform(0, 1, 400_000)                                   # true phase
    q = np.floor(u * 4).astype(int)
    f = u * 4 - q
    tap = np.searchsorted(edges, f, side="right") - 1                # 0..139 forward
    code = q * 256 + (10 + 139 - tap)                                # reversed numbering
    tab = tdc.code_density_table(code, reverse_taps=True)
    est = tab[code]
    centre = (q + (edges[tap] + edges[tap + 1]) / 2) / 4
    assert np.abs(est - centre).max() < 0.004
    tab_wrong = tdc.code_density_table(code, reverse_taps=False)
    assert np.abs(tab_wrong[code] - centre).max() > 0.05


# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("missing", [[], [700, 1500, 1501]])
def test_pulse_train(missing):
    P = 78.333
    n = 3000
    k = np.arange(n + len(missing))
    keep = np.setdiff1d(k, missing)
    rng = np.random.default_rng(3)
    S = 50 + P * keep + rng.normal(scale=0.003, size=keep.size)
    co = np.floor(S).astype(np.int64)
    pt = timing.build_pulse_train(co, S - co, divider=32, smooth_window=64)
    assert pt.n_missing_sync == len(missing)
    assert abs(pt.ticks_per_sync - P) < 1e-4
    # a time 10.25 pulses after sync period 1600 (after the gaps, if any)
    t = 50 + P * (1600 + 10.25 / 32)
    pulse, phase = pt.locate(np.array([t]))
    assert pulse[0] == 1600 * 32 + 10
    assert abs(phase[0] - 0.25) < 0.01


# --------------------------------------------------------------------------- #
def test_t3_roundtrip_with_ptufile(tmp_path):
    ptufile = pytest.importorskip("ptufile")
    rng = np.random.default_rng(5)
    n = 20000
    nsync = np.sort(rng.integers(0, 50_000_000, n))
    nsync[100:110] += 3_000_000_000 // 1000      # include a long gap
    nsync = np.sort(nsync)
    special = (rng.uniform(size=n) < 0.01).astype(np.uint8)
    channel = np.where(special == 1, rng.choice([1, 2, 4], n), rng.integers(0, 23, n)).astype(np.uint8)
    dtime = np.where(special == 1, 0, rng.integers(0, 5000, n)).astype(np.uint16)
    rec = encode_t3_fast(nsync, channel, dtime, special)
    assert np.array_equal(rec, encode_t3(nsync, channel, dtime, special))
    tags = [("File_GUID", "{00000000-0000-0000-0000-000000000000}"),
            ("Measurement_Mode", 3), ("Measurement_SubMode", 3),
            ("MeasDesc_GlobalResolution", 25e-9), ("MeasDesc_Resolution", 5e-12),
            ("TTResult_SyncRate", 40_000_000), ("TTResult_NumberOfRecords", 0),
            ("TTResultFormat_TTTRRecType", 0x01010304), ("TTResultFormat_BitsPerRecord", 32),
            ("ImgHdr_PixX", 10), ("ImgHdr_PixY", 10), ("ImgHdr_LineStart", 1),
            ("ImgHdr_LineStop", 2), ("ImgHdr_Frame", 3), ("ImgHdr_Dimensions", 3)]
    path = tmp_path / "t.ptu"
    write_ptu(str(path), tags, rec)
    p = ptufile.PtuFile(str(path))
    assert p.tags["TTResult_NumberOfRecords"] == rec.size
    d = p.decode_records()
    ev = (d["channel"] >= 0) | (d["marker"] > 0)
    d = d[ev]
    assert np.array_equal(d["time"].astype(np.int64), nsync)
    ph = special == 0
    assert np.array_equal(d["channel"][ph], channel[ph].astype(np.int8))
    assert np.array_equal(d["dtime"][ph], dtime[ph].astype(np.int16))
    assert np.array_equal(d["marker"][~ph], channel[~ph])


# --------------------------------------------------------------------------- #
@pytest.mark.skipif(not os.path.exists(os.path.join(ROOT, "inputs", "Matrix.raw")),
                    reason="input data not available")
def test_matrix_integration(tmp_path):
    """End-to-end on the real Matrix.raw with a temporary output directory."""
    import yaml
    from matrix2ptu.converter import convert
    from matrix2ptu.validate import validate
    cfg = yaml.safe_load(open(os.path.join(ROOT, "configs", "matrix.yaml")))
    cfg["input"]["raw"] = os.path.join(ROOT, cfg["input"]["raw"])
    cfg["input"]["msr"] = os.path.join(ROOT, cfg["input"]["msr"])
    cfg["output"]["directory"] = str(tmp_path)
    cp = tmp_path / "cfg.yaml"
    yaml.safe_dump(cfg, open(cp, "w"))
    res = convert(str(cp), keep_validation_data=True)
    out = validate(res, tif_path=os.path.join(ROOT, "inputs", "Intensity.tif"), run_flopa=True)
    failed = [c for c in out["checks"] if c["passed"] is False]
    assert not failed, failed


# --------------------------------------------------------------------------- #
# Block-wise processing must give exactly the same result as one block.
def test_decode_blockwise_identical(tmp_path, monkeypatch):
    path, _, _ = _synthetic_stream(tmp_path)
    whole = rawstream.decode(path)
    monkeypatch.setattr(rawstream, "BLOCK_WORDS", 97)
    blocks = rawstream.decode(path)
    for f in ("ph_element", "ph_code", "ph_coarse", "sy_code", "sy_coarse",
              "mk_bits", "mk_coarse", "mk_word"):
        assert np.array_equal(getattr(whole, f), getattr(blocks, f)), f
    assert whole.max_backstep == blocks.max_backstep


def test_pulse_train_blockwise_identical():
    rng = np.random.default_rng(7)
    keep = np.setdiff1d(np.arange(5000), [700, 1500, 1501, 4999 - 40])
    S = 50 + 78.333 * keep + rng.normal(scale=0.003, size=keep.size)
    a = timing.pulse_train_from_times(S.copy(), 32, 64, block=10 ** 9)
    b = timing.pulse_train_from_times(S.copy(), 32, 64, block=101)
    assert np.array_equal(a.sync_times, b.sync_times)
    assert np.array_equal(a.gap_pos, b.gap_pos) and np.array_equal(a.gap_cum, b.gap_cum)
    assert a.ticks_per_sync == b.ticks_per_sync and a.jitter_ticks_rms == b.jitter_ticks_rms


def test_write_ptu_t3_blockwise_identical(tmp_path):
    from matrix2ptu.ptu_writer import write_ptu_t3
    rng = np.random.default_rng(9)
    n = 5000
    nsync = np.sort(rng.integers(0, 5_000_000, n))
    nsync[3000:] += 2_000_000          # a gap needing more than 1023 wraps
    special = np.zeros(n, np.uint8)
    channel = rng.integers(0, 23, n).astype(np.uint8)
    dtime = rng.integers(0, 5000, n).astype(np.uint16)
    tags = [("File_GUID", "{0}"), ("TTResult_NumberOfRecords", 0)]
    write_ptu(str(tmp_path / "a.ptu"), tags, encode_t3_fast(nsync, channel, dtime, special))
    write_ptu_t3(str(tmp_path / "b.ptu"), tags, nsync, channel, dtime, special, block=333)
    assert (tmp_path / "a.ptu").read_bytes() == (tmp_path / "b.ptu").read_bytes()


def test_streamed_pulse_train_identical(tmp_path, monkeypatch):
    """Pulse train re-read block-wise from the file == in-memory version,
    with small blocks so that block edges, halos and gaps are exercised."""
    path, t_sync, _ = _synthetic_stream(tmp_path, n_sync=6000)
    monkeypatch.setattr(rawstream, "BLOCK_WORDS", 211)
    ev = rawstream.decode(path, keep_sync=True)
    # drop a few sync words to create gaps
    w = np.fromfile(path, dtype="<u4")
    is_sync = ((w >> 31) == 1) & (((w >> 26) & 0x1F) == 23)
    drop = np.flatnonzero(is_sync)[[700, 1500, 1501, 4000]]
    p2 = tmp_path / "gaps.raw"
    np.delete(w, drop).tofile(p2)
    ev = rawstream.decode(str(p2), keep_sync=True)
    rng = np.random.default_rng(2)
    table = np.sort(rng.uniform(0, 1, 1024))
    S = ev.sy_coarse.astype(np.float64) + table[ev.sy_code]
    ref = timing.pulse_train_from_times(S.copy(), 32, 64, block=97)
    st = timing.streamed_pulse_train(ev.sync, table, 32, 64, block=97)
    assert st.n_missing_sync == ref.n_missing_sync == 4
    assert st.ticks_per_sync == ref.ticks_per_sync and st.jitter_ticks_rms == ref.jitter_ticks_rms
    t = np.concatenate([rng.uniform(S[0] - 50, S[-1] + 50, 3000), ref.sync_times[::7]])
    for a, b in zip(ref.locate(t), st.locate(t)):
        assert np.array_equal(a, b)
    assert np.array_equal(ref.outside_mask(t), st.outside_mask(t))


def test_smoothed_block_chunking_identical():
    """The chunked moving mean must equal the one-shot version bit for bit."""
    rng = np.random.default_rng(4)
    n, half, block = 1000, 32, 300
    S = 50 + 78.333 * np.arange(n) + rng.normal(scale=0.003, size=n)
    raw = lambda a, b: S[a:b].copy()
    line = lambda a, b: 50 + 78.333 * np.arange(a, b)
    for j in range(4):
        one = timing._smoothed_block(raw, line, n, half, block, j, chunk=10 ** 9)
        for chunk in (1, 7, 64, 299):
            assert np.array_equal(one, timing._smoothed_block(raw, line, n, half, block, j, chunk=chunk))
    # and equal to the in-memory reference
    ref = S.copy()
    timing._smooth_inplace(ref, lambda a, b: S[a:b] - line(a, b), half, block)
    got = np.concatenate([timing._smoothed_block(raw, line, n, half, block, j, chunk=13) for j in range(4)])
    assert np.array_equal(ref, got)


def test_iter_ptu_records_blockwise(tmp_path):
    """Block-wise ptufile decoding with the overflow carry == one full decode,
    including blocks that end in overflow records and gaps > 1023 wraps."""
    ptufile = pytest.importorskip("ptufile")
    from matrix2ptu.validate import iter_ptu_records
    rng = np.random.default_rng(11)
    n = 3000
    nsync = np.sort(rng.integers(0, 3_000_000, n))
    nsync[2000:] += 5_000_000
    special = (rng.uniform(size=n) < 0.05).astype(np.uint8)
    channel = np.where(special == 1, 1, rng.integers(0, 23, n)).astype(np.uint8)
    dtime = np.where(special == 1, 0, rng.integers(0, 5000, n)).astype(np.uint16)
    tags = [("File_GUID", "{0}"), ("Measurement_Mode", 3), ("Measurement_SubMode", 3),
            ("MeasDesc_GlobalResolution", 25e-9), ("MeasDesc_Resolution", 5e-12),
            ("TTResult_SyncRate", 40_000_000), ("TTResult_NumberOfRecords", 0),
            ("TTResultFormat_TTTRRecType", 0x01010304), ("TTResultFormat_BitsPerRecord", 32)]
    path = str(tmp_path / "t.ptu")
    write_ptu(path, tags, encode_t3_fast(nsync, channel, dtime, special))
    p = ptufile.PtuFile(path)
    full = p.decode_records()
    p.close()
    for chunk in (5, 17, 1000, 10 ** 6):
        parts = list(iter_ptu_records(path, chunk))
        assert sum(d.size for _, d, _ in parts) == full.size
        for s0, d, t in parts:
            f = full[s0:s0 + d.size]
            ev = (f["channel"] >= 0) | (f["marker"] > 0)
            assert np.array_equal(t[ev], f["time"][ev].astype(np.int64))
            for name in ("dtime", "channel", "marker"):
                assert np.array_equal(d[name], f[name])


# --------------------------------------------------------------------------- #
# Streamed .obf input (event words inside an OBF container)
def _write_obf(path, stacks, file_version=2):
    """Minimal OBF writer for tests: stacks = [(name, array, footer_bytes, compression)]."""
    import struct
    import zlib
    from matrix2ptu.obf import FILE_MAGIC, STACK_MAGIC
    codes = {np.dtype(np.uint32): 0x10, np.dtype(np.int16): 0x08}
    out = bytearray(FILE_MAGIC + struct.pack("<IQ", file_version, 0) + b"\0" * 8)
    struct.pack_into("<Q", out, 14, len(out))
    for i, (name, arr, footer, comp) in enumerate(stacks):
        pos = len(out)
        data = arr.tobytes()
        if comp:
            data = zlib.compress(data)
        res = list(reversed(arr.shape)) + [1] * (15 - arr.ndim)
        h = STACK_MAGIC + struct.pack("<II", 1, arr.ndim) + struct.pack("<15I", *res)
        h += struct.pack("<15d", *([0.0] * 15)) + struct.pack("<15d", *([0.0] * 15))
        nb = name.encode()
        h += struct.pack("<5I", codes[arr.dtype], comp, 0, len(nb), 0)
        nxt_field = len(out) + len(h) + 16          # position of the 'next' field
        h += struct.pack("<QQQ", 0, len(data), 0) + nb
        out += h + data + footer
        if i + 1 < len(stacks):
            struct.pack_into("<Q", out, nxt_field, len(out))
    open(path, "wb").write(bytes(out))


def test_obf_stream_decodes_like_raw(tmp_path, monkeypatch):
    raw, _, _ = _synthetic_stream(tmp_path)
    words = np.fromfile(raw, dtype="<u4")           # including the zero padding
    obf = tmp_path / "syn.obf"
    _write_obf(obf, [("00000", words, b"<root><doc/></root>", 0)])
    from matrix2ptu.obf import find_event_stream
    es = find_event_stream(str(obf))
    assert es.n_words == words.size and es.stack_name == "00000"
    for block in (None, 97):
        if block:
            monkeypatch.setattr(rawstream, "BLOCK_WORDS", block)
        a, b = rawstream.decode(raw), rawstream.decode(str(obf))
        assert a.n_words == b.n_words and a.n_padding_words == b.n_padding_words == 7
        for f in ("ph_element", "ph_code", "ph_coarse", "sy_code", "sy_coarse",
                  "mk_bits", "mk_coarse", "mk_word"):
            assert np.array_equal(getattr(a, f), getattr(b, f)), f


def test_obf_without_event_stream_is_refused(tmp_path):
    from matrix2ptu.obf import find_event_stream
    img = tmp_path / "img.obf"
    _write_obf(img, [("STAR ORANGE multichannel", np.zeros((23, 5, 4), np.int16), b"", 1)])
    with pytest.raises(ValueError, match="no 1-D uint32 event-stream stack"):
        find_event_stream(str(img))
    comp = tmp_path / "comp.obf"
    _write_obf(comp, [("00000", np.arange(10, dtype=np.uint32), b"", 1)])
    with pytest.raises(ValueError, match="compressed"):
        find_event_stream(str(comp))


def test_tdc_gap_fill():
    tab = np.full(1024, np.nan)
    tab[10], tab[12], tab[300] = 0.1, 0.2, 0.3
    cal = tdc.TdcCalibration({0: tab.copy(), 23: tab.copy()}, "test")
    el = np.array([0, 0, 0, 0], np.uint8)
    code = np.array([10, 11, 11, 12], np.uint16)
    filled = cal.fill_gaps(el, code)
    assert filled == {0: 2}
    assert cal.tables[0][11] == pytest.approx(0.15)
    assert np.isnan(cal.tables[0][13])             # unused codes stay undefined


def test_config_needs_exactly_one_stream(tmp_path):
    from matrix2ptu.config import load_config
    for inp, ok in (({"raw": "a.raw"}, True), ({"obf": "a.obf"}, True),
                    ({"raw": "a.raw", "obf": "a.obf"}, False), ({"msr": "a.msr"}, False)):
        p = tmp_path / "c.yaml"
        p.write_text("input:\n" + "".join(f"  {k}: {v}\n" for k, v in inp.items()))
        if ok:
            load_config(str(p))
        else:
            with pytest.raises(ValueError, match="exactly one"):
                load_config(str(p))


@pytest.mark.skipif(not os.path.exists(os.path.join(ROOT, "inputs2", "MatrixOBF.obf")),
                    reason="input data not available")
def test_matrix_obf_integration(tmp_path):
    """End-to-end on the real streamed MatrixOBF.obf (metadata from the .obf itself,
    cross-checked against MatrixOBF.msr, which is also the validation reference)."""
    import yaml
    from matrix2ptu.converter import convert
    from matrix2ptu.validate import validate
    cfg = yaml.safe_load(open(os.path.join(ROOT, "configs", "matrix_obf.yaml")))
    for k in ("obf", "msr"):
        cfg["input"][k] = os.path.join(ROOT, cfg["input"][k])
    cfg["timing"]["tdc_calibration_file"] = os.path.join(ROOT, cfg["timing"]["tdc_calibration_file"])
    cfg["output"]["directory"] = str(tmp_path)
    cp = tmp_path / "cfg.yaml"
    yaml.safe_dump(cfg, open(cp, "w"))
    res = convert(str(cp), keep_validation_data=True)
    assert res["report"]["source"]["container"] == "obf"
    assert res["report"]["parameters"]["msr_vs_obf_parameters"]["value"] == "agree"
    out = validate(res, run_flopa=True)
    failed = [c for c in out["checks"] if c["passed"] is False]
    assert not failed, failed


# --------------------------------------------------------------------------- #
# pixel reassignment
def test_xcorr_shift_recovers_known_shift():
    from matrix2ptu import beadcal
    rng = np.random.default_rng(3)
    yy, xx = np.indices((128, 128))
    img = np.zeros((128, 128))
    for _ in range(25):                       # sparse field of Gaussian spots
        y0, x0 = rng.uniform(10, 118, 2)
        img += 100 * np.exp(-((yy - y0) ** 2 + (xx - x0) ** 2) / (2 * 3.0 ** 2))
    moved = beadcal.fourier_shift(img, 2.3, -1.7)
    sx, sy = beadcal.xcorr_shift(moved, img)
    assert abs(sx - 2.3) < 0.05 and abs(sy + 1.7) < 0.05


def test_gauss_fit_with_neighbour():
    from matrix2ptu import beadcal
    rng = np.random.default_rng(4)
    yy, xx = np.indices((25, 25), dtype=float)
    g = lambda y0, x0, a: a * np.exp(-((yy - y0) ** 2 + (xx - x0) ** 2) / (2 * 4.0 ** 2))
    img = rng.poisson(g(12.3, 11.6, 400) + g(12.0, 24.0, 400) + 2.0)
    p, err, ok = beadcal.fit_gauss(img, 12, 12, 4.0, neighbours=[(12.0, 24.0)])
    assert ok and abs(p[1] - 12.3) < 4 * err[1] + 0.05 and abs(p[2] - 11.6) < 4 * err[2] + 0.05
    p0, _, _ = beadcal.fit_gauss(img, 12, 12, 4.0)        # neighbour ignored -> pulled towards it
    assert p0[2] - 11.6 > p[2] - 11.6 + 0.2


def _split_ptu(path, photons, n_lines=4, acc=1, nx=10, dwell=10, mode="detector_resolved",
               seqs="0:CONF:x1"):
    """Minimal per-sequence PTU as the converter writes it: line windows of
    nx * dwell syncs; photons = [(element, line pass, x position in px, dtime)]."""
    from matrix2ptu.ptu_writer import write_ptu_t3
    ls = 1000 + np.arange(n_lines * acc) * 3 * nx * dwell
    le = ls + nx * dwell
    ev = [(0, 4, 0, 1)] + [(int(a), 1, 0, 1) for a in ls] + [(int(b), 2, 0, 1) for b in le]
    for e, k, x, d in photons:
        ev.append((int(ls[k] + x * dwell), e, d, 0))
    ev.sort(key=lambda r: (r[0], -r[3]))
    n, c, d, s = (np.array(v) for v in zip(*ev))
    tags = [("File_GUID", "{00000000-0000-0000-0000-000000000000}"), ("File_Comment", "test"),
            ("Measurement_Mode", 3), ("Measurement_SubMode", 3),
            ("MeasDesc_GlobalResolution", 25e-9), ("MeasDesc_Resolution", 5e-12),
            ("TTResult_SyncRate", 40_000_000), ("TTResult_NumberOfRecords", 0),
            ("TTResultFormat_TTTRRecType", 0x01010304), ("TTResultFormat_BitsPerRecord", 32),
            ("TTResultFormat_WrapAround", 1024),
            ("ImgHdr_PixX", nx), ("ImgHdr_PixY", n_lines * acc), ("ImgHdr_PixResol", 0.02),
            ("ImgHdr_LineStart", 1), ("ImgHdr_LineStop", 2), ("ImgHdr_Frame", 3), ("ImgHdr_Dimensions", 3),
            ("UsrMX_OutputMode", mode), ("UsrMX_Sequences", seqs)]
    write_ptu_t3(str(path), tags, n.astype(np.int64), c.astype(np.uint8), d.astype(np.uint16),
                 s.astype(np.uint8))


def test_reassign_ptu_moves_photons_keeps_timing(tmp_path):
    pytest.importorskip("ptufile")
    from matrix2ptu import ism
    from matrix2ptu.reassign import load_photons
    src = tmp_path / "src.ptu"
    # element 0 at (x 5.3, line 2), element 1 at (x 1.5, line 0), element 2 at (x 7.2, line 3)
    _split_ptu(src, [(0, 2, 5.3, 111), (1, 0, 1.5, 222), (2, 3, 7.2, 333)])
    shifts_px = np.zeros((23, 2))
    shifts_px[0] = (2.0, 1.2)       # -> x 3.3, line 2 - 1 = 1
    shifts_px[1] = (2.0, 0.0)       # -> x -0.5: leaves the field, dropped
    shifts_px[2] = (-1.5, -0.4)     # -> x 8.7, line 3 - 0 = 3
    cal = {"sequences": {"CONF": {"shift_nm": (shifts_px * 20).tolist()}}}
    rep = ism.reassign_ptu(str(src), cal, "CONF", str(tmp_path / "det.ptu"), str(tmp_path / "pool.ptu"))
    assert rep["stats"]["photons_dropped_out_of_field"] == 1
    out = load_photons(str(tmp_path / "det.ptu"), 1, 10, 4, 20)
    got = sorted((int(c), int(np.floor(x)), int(y), int(d)) for c, x, y, d in zip(out.ch, out.x, out.y, out.dtime))
    assert got == [(0, 3, 1, 111), (2, 8, 3, 333)]
    for c in ism.verify_reassigned(rep, run_tttrkit=False):
        assert c["structure"]["time_monotonic"] and c["structure"]["markers_same_as_source"]
        assert c["micro_time_multiset_per_channel_identical"] and c["image_vs_new_macro_times"]["identical"]
    pooled = load_photons(str(tmp_path / "pool.ptu"), 1, 10, 4, 20)
    assert set(pooled.ch.tolist()) == {0} and sorted(pooled.dtime.tolist()) == [111, 333]


def test_reassign_refuses_pooled_or_multi_sequence(tmp_path):
    pytest.importorskip("ptufile")
    from matrix2ptu import ism
    _split_ptu(tmp_path / "p.ptu", [(0, 0, 1, 1)], mode="pooled")
    with pytest.raises(ValueError, match="detector-resolved"):
        ism.scan_info(str(tmp_path / "p.ptu"))
    _split_ptu(tmp_path / "m.ptu", [(0, 0, 1, 1)], seqs="0:CONF:x1;1:STED:x3")
    with pytest.raises(ValueError, match="one PTU per sequence"):
        ism.scan_info(str(tmp_path / "m.ptu"))
