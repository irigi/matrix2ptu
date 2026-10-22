# MATRIX pixel reassignment: assessment of the available information

> **Update 2026-10-22.** Pixel reassignment with a bead calibration is now done; see [ISM_REPORT.md](ISM_REPORT.md). The shift scales here (nm/pitch) use each element versus the sum of the *other* elements; ISM_REPORT uses the pooled image as reference, which gives smaller numbers (e.g. 42 instead of 51 nm/pitch for the 640 nm data). The cross-correlation in `tools/analyze_array_geometry.py` used zero padding, which biases slightly towards zero; re-measured with the corrected estimator the values change by ≤ 2 %, so the conclusions here stand.

*2026-10-21. Assessment only: no pixel reassignment is implemented. Labels as in [REPORT.md](REPORT.md): [V] verified, [A] Abberior, [H] hypothesis, [U] unresolved.*

Analysis tool: `tools/analyze_array_geometry.py`.
- Results: `outputs/array_geometry.json` and `outputs/array_geometry_log.txt`.
- Figure: [f7_array_shifts.png](figures/f7_array_shifts.png).

## Short answer for the core facility

* **The element numbering and the array geometry are now established well enough to start [V].**
  * Data element *e* = Imspector label *e* + 1, in the 5-4-5-4-5 hexagonal layout of the Imspector screenshot.
  * The image-shift pattern of all seven datasets fits this layout. It does not fit any of 2000 random numberings: residual 4–16 nm vs ≥ 12–74 nm.
* **The Huygens values describe a lattice model, not a calibration.**
  * Huygens itself marks spacing, aspect ratio, rotation and mirroring as *"Derived value – verification needed"*.
  * They match our data in rotation (about 85°, with the opposite sign convention) and roughly in aspect ratio.
  * They do **not** contain the reassignment shifts: the scale factor between "spacing" and the actual shift is 0.5–0.65 in confocal and 0.15 in STED, and depends on the acquisition.
* **The Huygens CSV contains no shifts [V].** It is the photon total per element (23 numbers) of `MatrixLB.obf`. It is useful as an alignment and quality-control fingerprint only.
* **The Lightbox calibration cannot be transferred directly [V + H].**
  * `MatrixLB.obf` was recorded with an **8 AU** pinhole; the Imspector datasets with **1 AU**.
  * The measured shift scale differs by 25 % at the same wavelength and pixel size (63 vs 49 nm per element pitch).
  * The emission spot sits differently on the array.
* **Minimum to proceed:** shifts must be measured **from the data themselves** (adaptive pixel reassignment), or from a bead or dense-sample calibration recorded with **the same pinhole, wavelength and mode**. Both are possible with the current converter outputs, because the detector-resolved PTU keeps the element of every photon.

---

## 1. Files and what they contain

| file | content [V] | useful for reassignment? |
|---|---|---|
| `MatrixLayoutImspetor.PNG` | Imspector's Matrix1 panel: 23 circles labelled 1–23 in rows of 5, 4, 5, 4, 5 (hexagonal), "USB super speed" | **yes**: layout topology and label order (confirmed below). No dimensions, orientation in the sample, or optical scale. |
| `MatrixLB.obf` | Lightbox image OBF: one stack `STAR ORANGE multichannel`, 23 channels × 530 × 275 px, int16, 30 nm pixels, no time axis. OME-XML: EXC 561 nm, emission 616 nm (STAR ORANGE), UPLXAPO60XO NA 1.42 oil, pinhole **8.0 AU** (249 µm conjugate; Huygens: back-projected 2077 nm), dwell 50 µs, 2026-10-21 11:08:28 | **yes**: per-element images, so per-element shifts can be measured (§3). No element-position metadata in the file. |
| `MatrixLB_SVI_MP.PNG` | Huygens "Edit Microscopic Parameters" for MatrixLB: array detector model MATRIX/NanoSPAD, spacing 108.1 nm, aspect ratio (x:y) 0.7918, rotation 274.505°, mirroring none; **blue dots = "Derived value – verification needed"**; NA 1.42, embedding medium "ProLong Gold" 1.47 (orange/edited) | partly: a lattice model with an unknown convention (§4). |
| `MatrixLB_..._intensityTable_sum.csv` | Huygens "Array Detector Quality Control" export: header (model MATRIX/NanoSPAD, 23 detectors, 1 channel, background suppressed 0), one row `Channel 0` with 23 values. These equal the **per-element photon totals of `MatrixLB.obf`** (ratio 1.000 ± 0.003 at the CSV's 3 significant digits) | **no shifts**: intensity fingerprint only (§5). |
| detector-resolved PTUs of datasets A and B (this work) | every photon with its element 0–22, its time and its pixel | **yes**: the shifts can be measured on the data to be reassigned (§3). |

---

## 2. Geometry vs shifts: two different things

* **Geometric placement:** where element *e* sits on the chip (hexagonal pitch, orientation). This is a property of the detector and its mounting.
* **Reassignment shift:** how far the image recorded by element *e* is displaced in the sample, and therefore where its photons must be moved.
  * Roughly, shift ≈ factor × (element position back-projected to the sample).
  * The factor depends on the excitation and detection PSFs: wavelengths, Stokes shift, pinhole, and STED vs confocal (STED shrinks the shifts).
  * Real optics add aberrations, alignment offsets and pinhole clipping.
* So the geometry fixes the *pattern*, but the shifts must be measured, or modelled per acquisition setting.

---

## 3. Measured element shifts [V]

**Method:**
- For each dataset, 23 per-element images are formed. For the PTUs these use the image photons, per sequence.
- Each element image is cross-correlated with the sum of the *other* 22 elements (leave-one-out), with a sub-pixel peak fit.
- The 23 shift vectors are fitted with the Imspector layout using a general linear map (affine: scale, rotation, aspect, shear) and with a similarity map (scale, rotation only). Weights are photon counts.
- Rotation is given clockwise-positive in image coordinates (x to the right, y = scan line index downwards).

| dataset | photons | scale (nm per element pitch) | aspect (affine) | rotation | residual affine / similarity | best of 2000 random numberings |
|---|---|---|---|---|---|---|
| MatrixLB.obf (Lightbox, 561 nm, 8 AU) | 6.68 M | 62.6 | 0.74 | 83.4° | 14.6 / 19.9 nm | 74 nm |
| Orangebeads RAW, CONF (561 nm, 1 AU) | 4.97 M | 48.6 | 0.67 | 82.5° | 13.3 / 18.1 nm | 43 nm |
| Orangebeads OBF, CONF | 4.62 M | 49.4 | 0.64 | 82.6° | 13.0 / 18.8 nm | 43 nm |
| Orangebeads RAW, STED | 4.31 M | 13.4 | 0.57 | 83.1° | 4.6 / 6.5 nm | 12 nm |
| Orangebeads OBF, STED | 4.10 M | 14.2 | 0.57 | 82.4° | 4.4 / 6.6 nm | 13 nm |
| Matrix RAW, CONF (640 nm, 1 AU) | 0.54 M | 54.8 | 0.77 | 84.9° | 16.3 / 18.5 nm | 48 nm |
| Matrix OBF, CONF | 0.54 M | 51.1 | 0.75 | 85.2° | 16.2 / 18.5 nm | 45 nm |

**Findings:**
1. **Numbering [V]:** element = label − 1. In every dataset, none of 2000 random numberings fits as well as this one, by a factor of 3–5 in residual.
2. **Orientation is stable [V]:** 82–85° in all datasets, Lightbox and Imspector alike, RAW and OBF alike. No mirroring.
   - The 5-element rows of the layout run along the image **y** axis (scan-line direction).
   - Elements 0–4 lie at +x, elements 18–22 at −x.
3. **The shift scale depends on the acquisition [V]:**
   - confocal at 561 nm and 1 AU: 49 nm/pitch, reproducible between the RAW and OBF acquisitions to 2 %;
   - STED: 14 nm/pitch, i.e. 0.28 × confocal, as expected for a much smaller excitation PSF;
   - 640 nm confocal: 51–55 nm;
   - Lightbox at 8 AU: **63 nm**.
4. **Aspect ratio below 1 [V]:**
   - The affine fit is consistently better than a pure similarity fit, and the shift lattice is compressed along one axis.
   - The aspect ratio (0.57–0.77) differs between the datasets. A single detector "aspect ratio" is therefore not sufficient.
   - Possible contributions [H]: anisotropic PSFs (polarisation, STED pattern), scan-direction effects such as pixel dwell blur along x, and pinhole clipping.
5. **The outer elements deviate from the lattice [V]** (Fig. 7):
   - The weak edge elements (0, 3, 4, 18–22) are shifted *less* than a uniform lattice predicts, most strongly at 1 AU.
   - Likely reasons [H]: the 1 AU pinhole clips their light, and their low signal-to-background ratio pulls cross-correlation shifts towards zero.
   - A 4-parameter lattice is therefore not an accurate description of the shifts. Per-element measured shifts are needed.
6. **The emission spot is off-centre on the array at 561 nm [V]:**
   - Intensity-weighted centroid, in pitch units: (−0.38, −0.23) for orangebeads and (−0.27, −0.26) for Lightbox.
   - At 640 nm it is centred: (0.02, −0.01).
   - This indicates a wavelength-dependent alignment offset of about ⅓ of the element pitch [H: chromatic or alignment origin].
   - It makes the shift pattern asymmetric, so the reference point for reassignment has to be the measured centroid, not element 11.

---

## 4. The Huygens parameters

**Do they describe the same lattice?**
- The tool builds element positions from the Huygens numbers under 16 candidate conventions: rotation sign, axis of the aspect ratio, y flip, layout orientation. It then fits only a common scale factor *k* and an offset to the measured shifts.
- The best convention is the same in every dataset:
  - rotation applied with the **opposite sign** to our clockwise image convention, i.e. 274.505° ≡ −85.5°, a counter-clockwise rotation as displayed;
  - aspect ratio 0.7918 applied along the layout's 5-element rows;
  - no flip.
- With this convention the Huygens lattice explains the measured shifts **almost as well as a free affine fit**:
  - MatrixLB: 15.3 vs 14.6 nm residual;
  - orangebeads CONF: 14.9 vs 13.3 nm;
  - STED: 5.4 vs 4.6 nm.
- So rotation, aspect and mirroring are consistent with our measurements [V], under that inferred convention [H].

**What the "spacing" means [U]:**
- The best scale factors are *k* = 0.65 (MatrixLB), 0.51–0.52 (orangebeads CONF), 0.14–0.15 (STED) and 0.52–0.56 (640 nm).
- If 108.1 nm is the element pitch back-projected into the sample, the classic confocal ISM factor of about 0.5 fits the 1 AU data. It does **not** fit the Lightbox image that Huygens was given (0.65).
- That suggests Huygens took the spacing from a detector model and magnification, not from the measured shifts [H].
- With the Lightbox magnification to the pinhole plane (249 µm / 2.077 µm ≈ 120), 108.1 nm corresponds to about 13.0 µm on the detector side. This should be checked against the MATRIX element pitch and relay optics from Abberior or PI [U].

**Are the four parameters sufficient?**
- **For the element positions, nearly [V]:** with the inferred conventions, the hexagonal layout plus spacing, aspect, rotation and mirroring reproduces the shift *pattern* of the bright elements.
- **For the reassignment shifts, no:**
  - the scale factor *k* is not part of them and changes with pinhole, wavelength and STED;
  - the outer elements deviate;
  - the spot is off-centre;
  - the conventions (sign, origin, axis of the aspect ratio) are inferred, not documented [U].

---

## 5. The Huygens CSV

- **Rows and columns:** one row ("Channel 0", the only image channel); columns `Det. 0` … `Det. 22`, the detector elements in data order.
- **Values:** summed intensity per element over the whole image, in photon counts (no background suppression), 3 significant digits.
- **Comparison with other datasets:**
  - Lightbox (8 AU): the brightest elements are 11, 10 and 6 with 9–10 % each; the weakest is 22 with 0.6 %.
  - Imspector 1 AU data: they are more concentrated on the central elements (e.g. orangebeads element 22 has 0.06 %).
- **Use:** a quick alignment and pinhole fingerprint (where the spot sits on the array, how much the edge elements get). It is a reasonable QC record per session.
- **It contains no shifts and no positions.** If Huygens' Array QC tool can export *shifts* (the file name says `intensityTable_sum`, suggesting other tables exist), that export would be more useful [U].

---

## 6. Can MatrixLB's calibration be transferred?

**Not directly [V + H].**

| | MatrixLB.obf | dataset A | dataset B |
|---|---|---|---|
| software | Lightbox | Imspector | Imspector |
| excitation | 561 nm | 640 nm | 561 nm (+ STED775 in sequence 2) |
| pinhole | **8.0 AU** | **1.0 AU** (100 µm) | **1.0 AU** |
| pixel / dwell | 30 nm / 50 µs | 40 nm / 100 µs | 30 nm / 50 µs |
| measured shift scale | 63 nm/pitch | 51–55 | 49 (CONF), 14 (STED) |
| spot centroid (pitch units) | (−0.27, −0.26) | (0.02, −0.01) | (−0.38, −0.23) |

- **What transfers:** the numbering, the hexagonal topology, the orientation (about 83–85°) and no mirroring. These are the same in all datasets [V].
- **What does not transfer:** the shift magnitude, the outer-element behaviour and the spot position.
  - These depend on pinhole (8 vs 1 AU), wavelength, STED and alignment, and partly on the sample (signal-to-background ratio).
  - Pixel size changes only the units.
- **Optical alignment** can drift between sessions. The off-centre spot at 561 nm shows it matters.

---

## 7. What is needed to perform and validate reassignment

Prioritised, minimum first:

1. **Nothing new for a first implementation (already available):**
   - the detector-resolved PTU (or the `DecodedAcquisition` arrays) gives element, pixel and micro time for every photon;
   - **adaptive pixel reassignment** measures the 23 shifts on the image itself (as in §3) and applies them;
   - this works for bright, structured samples such as these datasets, per sequence (confocal and STED separately).
2. **A bead calibration per configuration:**
   - sub-resolution beads (e.g. 100 nm or smaller, matching dyes), recorded with **the same pinhole (1 AU or the setting used), wavelengths, objective and mode (CONF and STED)** as the experiments, at about 20–30 nm pixels;
   - saved as `.raw` or streamed `.obf` **plus** `.msr`, ideally on the same day as the samples;
   - this gives well-determined shifts for weak edge elements (the main weakness of data-derived shifts), and a ground truth for validation.
3. **Documentation from Abberior / PI Imaging:**
   - physical element pitch and layout orientation;
   - the optical magnification to the detector;
   - whether the MATRIX is placed behind the pinhole, and how the pinhole setting affects the outer elements;
   - whether Imspector/Lightbox store element positions anywhere.
4. **From Huygens (if available):**
   - the per-element shift table of its Array QC tool, and its definitions of spacing, aspect ratio and rotation;
   - Huygens' reassigned image of `MatrixLB.obf` and of a converted dataset, to serve as an external reference.
5. **A pinhole series on beads** (e.g. 1 AU and 8 AU at the same position): this would confirm or refute the pinhole hypothesis [H] for the difference between Lightbox and Imspector shifts.

**Validation of a reassignment** (planned, not done):
- beads: FWHM before and after vs the expected √2 gain (confocal), and FRC resolution on structured samples;
- closure tests:
  - shifts from the two separate acquisitions of the same sample agree (here 49 vs 49 nm/pitch, RAW vs OBF);
  - the STED/CONF shift ratio is consistent;
  - reassigned element images correlate better with each other than unshifted ones;
- the photon count is conserved;
- the lifetime per pixel is unchanged on average (reassignment moves photons, not their micro times; a check for errors).

## 8. What can be implemented now, and what needs more

| step | available now? |
|---|---|
| per-element images and shifts from converted data (confocal and STED) | **yes** (§3 tool) |
| adaptive pixel reassignment with data-derived shifts, intensity and FLIM (photon-level, micro times kept) | **yes**, implementable on `DecodedAcquisition` / detector-resolved PTU |
| a lattice model with Huygens-compatible parameters for comparison | **yes** (conventions inferred, not documented) |
| reliable shifts for the dim edge elements | needs **bead calibration** per configuration |
| absolute geometric calibration (physical pitch, magnification) | needs **Abberior/PI documentation** |
| a transfer rule between pinhole settings or wavelengths | needs **bead series** (pinhole, wavelength) |
| independent validation | needs beads, and ideally the Huygens shift export and reassigned references |
