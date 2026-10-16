# matrix2ptu — Abberior MATRIX `.raw` → PicoQuant `.ptu`

Converts the photon stream that Imspector writes for the Abberior **MATRIX**
detector (a 23-element PI Imaging SPAD23 array with per-element TDCs) into
**event-based PicoQuant PTU files** (one record per photon, with macro time,
micro time and detector channel, plus scan markers).

* **detector-resolved** PTU: channel = detector element 0…22
* **pooled** PTU: all elements in channel 0 (pooled per photon, after per-element delay correction)
* optional: one PTU per acquisition sequence (for software that cannot handle line-interleaved sequences)

Read **[docs/REPORT.md](docs/REPORT.md)** for what was established about the
format, how the conversion works, the validation results and the limitations.

---

## Step-by-step guide

1. **Install once** (details in §1):
   ```bash
   python3 -m venv .venv && source .venv/bin/activate
   pip install -e ".[validate,dev]"
   ```
2. **Collect the input files** of one measurement: the `.raw` file (required)
   and the `.msr` file of the same measurement (strongly recommended). The
   `.raw` has no header: without the `.msr`, image size, dwell time, sequence
   structure, laser frequency and per-element delays must be written in the
   config.
3. **Inspect the raw file** (optional, helps to fill in the config):
   ```bash
   matrix2ptu inspect inputs/Matrix.raw --msr inputs/Matrix.msr
   ```
4. **Write a YAML config for this dataset.** Copy the closest example,
   [`configs/matrix.yaml`](configs/matrix.yaml) (one sequence) or
   [`configs/sequence.yaml`](configs/sequence.yaml) (line-interleaved
   sequences), and edit at least:
   * `input.raw`, `input.msr` and `output.directory`
   * `output.modes`, and `output.split_sequences: true` if you need one PTU per sequence (e.g. for SymPhoTime)
   * `scan.*`: image size, pixel size, dwell **per line pass**
   * `acquisition.sequences`: name and `line_accumulation` of each sequence, in acquisition order
   * leave the `timing.*` delays and laser frequency as `auto` when the `.msr` is available

   All keys are explained in §3.
5. **Dry run:** checks all parameters and writes only the report, no PTU:
   ```bash
   matrix2ptu convert configs/mydata.yaml --dry-run
   ```
   Read the printed values, their sources and any warnings, and fix the config if needed.
6. **Convert and validate:**
   ```bash
   matrix2ptu convert configs/mydata.yaml --validate
   # single-sequence data with a reference intensity image:
   matrix2ptu convert configs/matrix.yaml --validate --tif inputs/Intensity.tif
   ```
   `--no-flopa` skips the tttrkit/FLOPA check. Exit code 0 = ok, 1 = error,
   2 = a validation check failed.
7. **Check the outputs** in `output.directory` (table in §2), especially
   `*_conversion_report.json` and `validation_report.json`.
8. **Open the PTU files** in FLOPA or SymPhoTime (§4). In FLOPA set the line
   accumulations to match the config, e.g. `(1)` for Matrix, `(2, 5)` for Sequence.

---

## 1. Installation (Linux/macOS/Windows, Python ≥ 3.10)

```bash
cd MaruskaConvertor
python3 -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
pip install -e ".[validate,dev]"   # core needs only numpy + pyyaml
```

`[validate]` installs the independent readers used for checking the output
(`ptufile`, `tttrkit` = FLOPA's reader, `tifffile`, `xarray`, `scipy`, `matplotlib`).

## 2. Running

```bash
# look at a raw file (and the .msr of the same measurement)
matrix2ptu inspect inputs/Matrix.raw --msr inputs/Matrix.msr

# convert, then validate against the source data and references
matrix2ptu convert configs/matrix.yaml   --validate --tif inputs/Intensity.tif
matrix2ptu convert configs/sequence.yaml --validate

# check settings without writing PTU files
matrix2ptu convert configs/sequence.yaml --dry-run
```

(Without installing: `PYTHONPATH=src python -m matrix2ptu ...`.)

Outputs go to the directory given in the config, e.g. `outputs/Matrix/`:

| file | content |
|---|---|
| `Matrix_detector_resolved.ptu` | 23 detector channels |
| `Matrix_pooled.ptu` | 1 channel |
| `Matrix_conversion_report.json` | every parameter used **and where it came from** (config / .msr / measured), warnings, counts |
| `validation_report.json` | results of all checks (with `--validate`) |
| `Matrix_tdc_calibration.npz` | the TDC calibration tables used (re-usable via `timing.tdc_calibration_file`) |

Exit code: 0 = ok, 1 = error (message explains what to fix), 2 = a validation check failed.
Set `MATRIX2PTU_DEBUG=1` to see full Python tracebacks.

Runtime/memory (this PC): Matrix.raw (90 MB) 2 s / 1.3 GB; Sequence.raw (600 MB) 20 s / 7 GB.
Memory scales with file size (≈12× the .raw size); `--validate` adds about 1 GB and ~50 s for
Sequence.raw. Planned reductions: [docs/BACKLOG.md](docs/BACKLOG.md).

## 3. Configuration

One YAML file per dataset — settings for one acquisition are never applied
to another. Commented examples: [`configs/matrix.yaml`](configs/matrix.yaml),
[`configs/sequence.yaml`](configs/sequence.yaml).

Every value can be written explicitly (used as given, then cross-checked) or
set to `auto` (taken from the `.msr` file of the same measurement or measured
from the raw stream). The program prints the final value and its source for
every parameter, and warns about disagreements.

| key | meaning | can it be inferred? |
|---|---|---|
| `scan.pixels_x`, `pixels_y` | image size | `.msr` (checked against the number of line markers) |
| `scan.pixel_size_um` | pixel size (header only) | `.msr` |
| `scan.pixel_dwell_us` | dwell **per line pass** | `.msr` (`pixel_steps/total`); checked against the line period |
| `scan.dwell_convention` | only `per_line_pass` (the dwell is never multiplied by the accumulation) | — |
| `scan.bidirectional`, `scan_type` | only monodirectional linear scanning is implemented | `.msr` |
| `scan.frames` | number of frames | checked against frame markers |
| `scan.line_start_offset_us` | shift of the pixel grid relative to the line marker | **user** (default 0) |
| `acquisition.sequences` | ordered list of `{name, line_accumulation}` — the order in which they are acquired for each y-line | `.msr` line steps (checked against the raw photon counts per pass) |
| `detector.n_elements`, `sync_element`, `sync_divider` | 23, 23, 32 (established from the data) | — |
| `timing.laser_frequency_hz` | laser repetition rate | measured: pulses per line × line frequency (`.msr`) — **needs the .msr or a user value** |
| `timing.microtime_resolution_ps` | PTU micro-time bin (default 5 ps) | user choice |
| `timing.element_delays_ps`, `common_delay_ps` | per-element timing offsets | `.msr` (Imspector calibration) — **needed for the pooled output** |
| `output.modes` | `detector_resolved`, `pooled` | — |
| `output.split_sequences` | also write one PTU per sequence | — |
| `output.hw_type` | PTU `HW_Type` string | — |

**Settings that cannot be inferred from the `.raw` alone** (the file has no
header): image size, dwell time, sequence/accumulation structure, laser
frequency and per-element delays. Supply the `.msr` of the same measurement or
write them in the config.

## 4. Opening the PTU files

**FLOPA** (napari-flopa): open the PTU, set *lines/pixels* and
**line accumulations** — Matrix: `(1)`; Sequence: `(2, 5)` (CONF then STED).
Markers: frame = 4, line start = 1, line stop = 2 (FLOPA defaults).
The detector-resolved file gives a `channel` dimension of 23.

**SymPhoTime 64**: not tested yet (see REPORT §8). Start with
`Matrix_pooled.ptu`. For Sequence data use `split_sequences: true`; each
per-sequence file then has `pixels_y × line_accumulation` lines (consecutive
passes of the same y-line), because SymPhoTime cannot combine accumulations.

## 5. Code layout (for later extensions such as pixel reassignment)

```
src/matrix2ptu/
  rawstream.py   decode 32-bit MATRIX words, unwrap coarse counter
  tdc.py         TDC code-density calibration
  timing.py      laser pulse train from the /32 sync; macro/micro time
  msr_meta.py    read acquisition metadata from .msr/.obf  (obf.py = container reader)
  config.py      YAML settings, provenance, cross-checks
  converter.py   pipeline: DecodedAcquisition (all per-photon quantities) + PTU assembly
  ptu_writer.py  PTU header + HydraHarp-v2 T3 record encoder
  validate.py    checks with independent readers (ptufile, tttrkit)
tools/           figures, FLOPA headless check
tests/           unit tests (synthetic) + integration test on Matrix.raw
```

`converter.DecodedAcquisition` holds, per photon: detector element, absolute
time, laser-pulse index, micro time, line pass and sequence — the natural
starting point for pixel reassignment (which needs the element geometry; see
REPORT §9).

Run tests: `pytest -q`.
