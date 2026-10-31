# matrix2ptu

**Convert Abberior MATRIX photon streams (`.raw` or streamed `.obf`) into
event-based PicoQuant `.ptu` files**, ready for FLIM analysis in FLOPA,
SymPhoTime or any PTU reader. Optionally, apply bead-calibrated pixel
reassignment (ISM) to the converted data.

The Abberior **MATRIX** detector is a 23-element SPAD array (PI Imaging SPAD23)
with a time-to-digital converter (TDC) per element. Imspector can stream its
photons to disk, but only in its own formats. `matrix2ptu` decodes that stream
and writes one PTU record per photon (macro time, micro time, detector channel)
plus frame and line markers.

**Contents:**
[Features](#features) ·
[Installation](#installation) ·
[Quick start](#quick-start) ·
[Commands](#commands) ·
[Configuration](#configuration) ·
[Output files](#output-files) ·
[Opening the PTU files](#opening-the-ptu-files) ·
[Pixel reassignment](#pixel-reassignment-ism) ·
[Performance](#performance) ·
[Documentation](#documentation) ·
[Development](#development)

---

## Features

* **Two input formats**
  * headerless **`.raw`** stream (Imspector Matrix1 streaming format "raw"),
    with the acquisition parameters taken from the `.msr` of the same measurement;
  * **streamed `.obf`** (Matrix1 streaming format "OBF"), which holds the same
    event words *and* the acquisition parameters.
* **Two kinds of output**
  * **detector-resolved** PTU: channel = detector element 0…22;
  * **pooled** PTU: all elements in channel 0, after per-element delay correction.
* **Line-interleaved sequences** (e.g. CONF + STED with different line
  accumulations), optionally split into one PTU per sequence for software that
  cannot handle interleaved sequences.
* **TDC calibration** (code-density) per element, re-usable for dim data.
* **Provenance for every parameter:** each value is reported together with where
  it came from (config, `.msr`/`.obf`, or measured from the stream), and
  disagreements are flagged as warnings.
* **Validation** of the output with independent readers (`ptufile`, FLOPA's
  `tttrkit`) against the source data and Imspector's own images.
* **Pixel reassignment (ISM)** with per-element shifts calibrated on beads.

Image-type files (`.msr`, Lightbox image `.obf`) contain no photon events and
cannot be converted to event-based PTU; see [docs/OBF_REPORT.md](docs/OBF_REPORT.md).

---

## Installation

Requires **Python 3.12 or newer** (the converter core runs on 3.10+, but the
validation readers need 3.12). The core depends only on `numpy` and `pyyaml`.

### Windows 11

1. **Install Python** from <https://www.python.org/downloads/windows/>
   (the "Windows installer (64-bit)"). On the first installer screen, tick
   **"Add python.exe to PATH"**.
2. **Get matrix2ptu**, either
   * on the GitHub page of this repository: **Code → Download ZIP**, then
     unzip it to a folder of your choice (for example `C:\matrix2ptu`), or
   * with Git: `git clone https://github.com/irigi/matrix2ptu.git`
3. **Double-click `install_windows.bat`** in that folder. It creates a private
   Python environment (`.venv`) inside the folder and installs everything into
   it; nothing else on the computer is changed. It ends with
   *"matrix2ptu is installed."*
4. **Double-click `matrix2ptu_terminal.bat`**. It opens a terminal in the
   folder where the `matrix2ptu` command works:
   ```bat
   matrix2ptu --help
   matrix2ptu convert configs\my_data.yaml --dry-run
   ```

Updating: replace the folder contents with a new download (or `git pull`) and
run `install_windows.bat` again. If you move the folder, run
`install_windows.bat -Recreate` from a terminal in the new place.

Installer options (from a terminal: `install_windows.bat -Dev`, etc.):

| option | effect |
|---|---|
| `-Dev` | also install `pytest`, to run the tests |
| `-NoValidate` | install only the converter, without the validation readers |
| `-Recreate` | delete `.venv` and build it again from scratch |

Paths in config files may use either `/` or `\`. If a path contains spaces or
starts with a drive letter, writing it in single quotes is safest:
`raw: 'D:\Data 2026\Matrix.raw'` (do not use double quotes around Windows paths;
YAML treats `\` inside double quotes as an escape character).

### Linux / macOS

```bash
git clone https://github.com/irigi/matrix2ptu.git
cd matrix2ptu
python3 -m venv .venv
source .venv/bin/activate
pip install -e ".[validate,dev]"     # or: pip install -e .   (converter only)
```

`[validate]` installs the independent readers used for checking the output
(`ptufile`, `tttrkit` = FLOPA's reader, `tifffile`, `xarray`, `scipy`,
`matplotlib`); `[dev]` installs `pytest`. The exact versions used for the
reports are pinned in [requirements.txt](requirements.txt).

Without installing, the package also runs as `PYTHONPATH=src python -m matrix2ptu ...`.

---

## Quick start

1. **Collect the input files** of one measurement:
   * the event stream (required): the `.raw` file **or** the streamed `.obf` file;
   * the `.msr` file of the same measurement (strongly recommended).

   A `.raw` file has no header. Without the `.msr`, image size, dwell time,
   sequence structure, laser frequency and per-element delays must be written
   in the config. A streamed `.obf` carries these parameters itself; its `.msr`
   is then used for cross-checks and as the validation reference.

2. **Inspect the stream** (optional, helps to fill in the config):
   ```bash
   matrix2ptu inspect data/Matrix.raw --msr data/Matrix.msr
   matrix2ptu inspect data/MatrixOBF.obf --msr data/MatrixOBF.msr
   ```

3. **Write a YAML config for this dataset.** Copy the closest example from
   [`configs/`](configs/) and edit it:

   | example | use for |
   |---|---|
   | [`matrix.yaml`](configs/matrix.yaml) | `.raw`, one sequence |
   | [`sequence.yaml`](configs/sequence.yaml) | `.raw`, line-interleaved sequences |
   | [`matrix_obf.yaml`](configs/matrix_obf.yaml) | streamed `.obf`, one sequence |
   | [`orangebeads_obf.yaml`](configs/orangebeads_obf.yaml) | streamed `.obf`, CONF + STED |
   | [`beads_obf_split.yaml`](configs/beads_obf_split.yaml) | streamed `.obf`, one PTU per sequence (input for pixel reassignment) |

   Edit at least:
   * `input.raw` **or** `input.obf`, then `input.msr` and `output.directory`;
   * `scan.*`: image size, pixel size, dwell **per line pass**;
   * `acquisition.sequences`: name and `line_accumulation` of each sequence, in
     acquisition order;
   * `output.modes`, and `output.split_sequences: true` if you need one PTU per
     sequence (e.g. for SymPhoTime);
   * for dim data (fewer than about 1000 photons in some element): set
     `timing.tdc_calibration_file` to a table from a bright acquisition of the
     same session (see [OBF_REPORT §4.2](docs/OBF_REPORT.md));
   * leave the `timing.*` delays and the laser frequency as `auto` when the
     `.msr` is available.

   The example configs point to the development datasets (`inputs/`,
   `inputs2/`, …), which are not part of the repository. Relative paths are
   relative to the folder you run `matrix2ptu` from.

4. **Dry run**: checks all parameters and writes only the report, no PTU:
   ```bash
   matrix2ptu convert configs/my_data.yaml --dry-run
   ```
   Read the printed values, their sources and any warnings, and fix the config
   if needed.

5. **Convert and validate:**
   ```bash
   matrix2ptu convert configs/my_data.yaml --validate
   # single-sequence data with a reference intensity image from Imspector:
   matrix2ptu convert configs/my_data.yaml --validate --tif data/Intensity.tif
   ```

6. **Check the outputs** in `output.directory` ([table below](#output-files)),
   especially `*_conversion_report.json` and `validation_report.json`.

7. **Open the PTU files** in FLOPA or SymPhoTime ([below](#opening-the-ptu-files)).

---

## Commands

| command | purpose |
|---|---|
| `matrix2ptu inspect STREAM [--msr MSR]` | summarise a `.raw` / streamed `.obf` (and its metadata) |
| `matrix2ptu convert CONFIG [--dry-run] [--validate] [--tif TIF] [--no-flopa]` | convert according to a YAML config |
| `matrix2ptu ism-calibrate --seq NAME=PTU ... --detect-in NAME --out CAL.json [--metadata OBF]` | per-element shifts from a bead acquisition |
| `matrix2ptu ism-reassign PTU --calibration CAL.json --sequence NAME --out-dir DIR [--metadata OBF]` | apply a shift calibration to one sequence |

* `--dry-run`: decode and check the parameters, write only the report.
* `--validate`: run the validation checks after the conversion.
* `--no-flopa`: skip the tttrkit/FLOPA reconstruction check.

**Exit codes:** 0 = ok, 1 = error (the message explains what to fix),
2 = a validation check failed. Set the environment variable
`MATRIX2PTU_DEBUG=1` to see full Python tracebacks.

---

## Configuration

One YAML file per dataset; settings for one acquisition are never applied to
another. Every value can be written explicitly (used as given, then
cross-checked) or set to `auto` (taken from the `.msr`/`.obf` of the same
measurement, or measured from the stream). The program prints the final value
and its source for every parameter.

```yaml
input:
  raw: data/Matrix.raw            # or  obf: data/MatrixOBF.obf
  msr: data/Matrix.msr            # optional, strongly recommended
output:
  directory: outputs/Matrix
  modes: [detector_resolved, pooled]
  split_sequences: false
scan:
  pixels_x: 400
  pixels_y: 400
  pixel_size_um: 0.05
  pixel_dwell_us: 100             # per line pass
acquisition:
  sequences:
    - name: Alexa488_CONF
      line_accumulation: 1
timing:
  laser_frequency_hz: auto
  element_delays_ps: auto
  common_delay_ps: auto
```

| key | meaning | can it be inferred? |
|---|---|---|
| `input.raw` / `input.obf` | the event stream: a `.raw` file or a streamed `.obf` (exactly one) | — |
| `input.msr` | metadata and validation reference | optional for `.obf` (its own parameters are used) |
| `scan.pixels_x`, `pixels_y` | image size | `.msr` (checked against the number of line markers) |
| `scan.pixel_size_um` | pixel size (header only) | `.msr` |
| `scan.pixel_dwell_us` | dwell **per line pass** | `.msr` (`pixel_steps/total`); checked against the line period |
| `scan.dwell_convention` | only `per_line_pass` (the dwell is never multiplied by the accumulation) | — |
| `scan.bidirectional`, `scan_type` | only monodirectional linear scanning is implemented | `.msr` |
| `scan.frames` | number of frames | checked against frame markers |
| `scan.line_start_offset_us` | shift of the pixel grid relative to the line marker | **user** (default 0) |
| `acquisition.sequences` | ordered list of `{name, line_accumulation}`, in the order they are acquired for each y-line | `.msr` line steps (checked against the photon counts per pass) |
| `detector.n_elements`, `sync_element`, `sync_divider` | 23, 23, 32 (established from the data) | — |
| `timing.laser_frequency_hz` | laser repetition rate | measured: pulses per line × line frequency — **needs the `.msr` or a user value** |
| `timing.microtime_resolution_ps` | PTU micro-time bin (default 5 ps) | user choice |
| `timing.element_delays_ps`, `common_delay_ps` | per-element timing offsets | `.msr` (Imspector calibration) — **needed for the pooled output** |
| `timing.tdc_calibration_file` | re-use a TDC table (`*_tdc_calibration.npz`) | needed for dim data; missing TDC codes are filled from neighbouring codes |
| `output.modes` | `detector_resolved`, `pooled` | — |
| `output.split_sequences` | also write one PTU per sequence | — |
| `output.hw_type` | PTU `HW_Type` string | — |

Fully commented examples: [`configs/matrix.yaml`](configs/matrix.yaml),
[`configs/sequence.yaml`](configs/sequence.yaml).

---

## Output files

Written to `output.directory`, e.g. `outputs/Matrix/`:

| file | content |
|---|---|
| `Matrix_detector_resolved.ptu` | 23 detector channels |
| `Matrix_pooled.ptu` | 1 channel |
| `Matrix_conversion_report.json` | every parameter used **and where it came from** (config / .msr / measured), warnings, counts |
| `validation_report.json` | results of all checks (with `--validate`) |
| `Matrix_tdc_calibration.npz` | the TDC calibration tables used (re-usable via `timing.tdc_calibration_file`) |

With `split_sequences: true`, one PTU per sequence is written in addition.

---

## Opening the PTU files

**FLOPA** (napari-flopa): open the PTU and set *lines/pixels* and the
**line accumulations** to match the config, e.g. `(1)` for a single sequence or
`(2, 5)` for CONF with 2 and STED with 5 accumulations. Markers: frame = 4,
line start = 1, line stop = 2 (FLOPA defaults). The detector-resolved file gives
a `channel` dimension of 23. Photons recorded during the scanner flyback stay in
the PTU but fall outside the line windows; image readers ignore them, whole-file
decay histograms include them ([OBF_REPORT §6.3](docs/OBF_REPORT.md)).

**SymPhoTime 64**: not tested yet ([REPORT §8](docs/REPORT.md)). Start with the
pooled PTU. For sequence data use `split_sequences: true`; each per-sequence
file then has `pixels_y × line_accumulation` lines (consecutive passes of the
same y-line), because SymPhoTime cannot combine accumulations.

---

## Pixel reassignment (ISM)

The 23 detector elements see the sample from slightly shifted positions.
Pixel reassignment moves each photon by its element's shift, calibrated on an
acquisition of isolated beads with the same settings
([docs/ISM_REPORT.md](docs/ISM_REPORT.md)).

```bash
# 1. convert beads and target with one PTU per sequence (output.split_sequences: true)
matrix2ptu convert configs/beads_obf_split.yaml
matrix2ptu convert configs/dense_obf_split.yaml

# 2. per-element shifts from isolated beads, one table per sequence
matrix2ptu ism-calibrate --seq CONF=<beads seq0 detector-resolved PTU> \
    --seq STED=<beads seq1 detector-resolved PTU> --detect-in STED \
    --metadata data/Beads.obf --out my_calibration.json

# 3. apply to one sequence of the target
#    (writes detector-resolved + pooled reassigned PTUs, checks and a report)
matrix2ptu ism-reassign <target seq0 detector-resolved PTU> \
    --calibration my_calibration.json --sequence CONF \
    --out-dir outputs/reassigned --metadata data/Dense.obf
```

(On Windows, write each command on one line, or end lines with `^` instead of `\`.)

Reassigned PTUs keep every photon's micro time and detector channel. The macro
time of a moved photon encodes its new position, so it is no longer the
detection time; use the original PTUs for anything that needs photon detection
times. An example calibration is in
[`configs/calibration/`](configs/calibration/).

---

## Performance

Measured on a desktop PC (Linux): Matrix.raw (90 MB) 3 s / 0.44 GB peak memory;
Sequence.raw (600 MB) 31 s / 1.5 GB. Conversion needs about 2.5× the `.raw`
size in memory, mostly the photons (about 45 bytes per photon at the peak).
With `--validate`: 5 s / 1.2 GB and 52 s / 2.7 GB. See
[docs/BACKLOG.md](docs/BACKLOG.md) for planned improvements.

---

## Documentation

| report | content |
|---|---|
| [docs/REPORT.md](docs/REPORT.md) | the `.raw` stream format, how the conversion works, validation, limitations |
| [docs/OBF_REPORT.md](docs/OBF_REPORT.md) | streamed `.obf` files, `.obf` → PTU, validation on 640 nm and 561 nm (CONF + STED) data, RAW vs OBF |
| [docs/PIXEL_REASSIGNMENT.md](docs/PIXEL_REASSIGNMENT.md) | what the available geometry and shift information supports for pixel reassignment |
| [docs/ISM_REPORT.md](docs/ISM_REPORT.md) | bead calibration and pixel reassignment: method, validation, results |
| [docs/BACKLOG.md](docs/BACKLOG.md) | open items and ideas |

---

## Development

```
src/matrix2ptu/
  rawstream.py   decode 32-bit MATRIX words, unwrap the coarse counter
  tdc.py         TDC code-density calibration
  timing.py      laser pulse train from the /32 sync; macro/micro time
  msr_meta.py    acquisition metadata from .msr/.obf
  obf.py         .obf container reader, incl. locating the event stack of a streamed .obf
  config.py      YAML settings, provenance, cross-checks
  converter.py   pipeline: DecodedAcquisition (all per-photon quantities) + PTU assembly
  ptu_writer.py  PTU header + HydraHarp-v2 T3 record encoder
  validate.py    checks with independent readers (ptufile, tttrkit)
  beadcal.py     bead detection, per-element Gaussian fits, shift calibration
  reassign.py    photon positions from a per-sequence PTU; event-level reassignment
  ism.py         calibration / reassignment workflow, calibration files, checks
tools/           figures, FLOPA headless and FLIM checks, RAW-vs-OBF comparison,
                 .obf container check, detector-geometry analysis, ISM round scripts
tests/           unit tests (synthetic data) + integration tests on the datasets
```

`converter.DecodedAcquisition` holds, per photon (sorted by laser pulse): detector
element, laser-pulse index, micro time and line pass (sequence via `seq_of_pass`);
with `keep_validation_data=True` also the absolute time and the phase. Pixel
reassignment works on the per-sequence detector-resolved PTUs instead
(`reassign.load_photons`), so it can be applied to any converted file.

**Tests:** `pytest -q`. Tests that need the measurement datasets are skipped
when the data are not present. GitHub Actions runs the tests on Linux and on
Windows (the Windows job installs with `install_windows.ps1`).

Additional checks used for the reports ([OBF_REPORT §5–6](docs/OBF_REPORT.md)):

```bash
# the .obf container route gives byte-identical PTU records to the same words as .raw
python tools/obf_container_check.py configs/matrix_obf.yaml
# FLOPA's per-pixel counts, mean arrival times and phasors vs Imspector's cube
python tools/flopa_flim_check.py outputs/MatrixOBF/MatrixOBF_pooled.ptu inputs2/MatrixOBF.msr 125 125 1
# compare two acquisitions of the same area (e.g. RAW vs OBF)
python tools/compare_pair.py outputs/MatrixRAW/MatrixRAW outputs/MatrixOBF/MatrixOBF
# per-element image shifts / detector geometry
python tools/analyze_array_geometry.py
# the whole pixel-reassignment round (conversion, calibration, transfer test, reassignment, figures)
python tools/ism_inputs3.py all
```

The FLOPA headless check additionally needs `napari-flopa` (see
[requirements.txt](requirements.txt)).

---

Author: Jan Olšina. All rights reserved (no license has been chosen yet).
