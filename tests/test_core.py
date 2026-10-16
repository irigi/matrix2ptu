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
    res = convert(str(cp))
    out = validate(res, tif_path=os.path.join(ROOT, "inputs", "Intensity.tif"), run_flopa=True)
    failed = [c for c in out["checks"] if c["passed"] is False]
    assert not failed, failed
