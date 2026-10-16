# Backlog

Planned work that is agreed but not scheduled. See `REPORT.md` §9 for open questions.

## Done

### A. Reduce memory use: smaller and fewer arrays  *(done 2026-10-16)*

**Problem.** Peak memory was about 12× the `.raw` size. Most of it was not photons:
sync records are 91 % (Matrix) and 77 % (Sequence) of all words, and several
full-length 8-byte arrays were kept at the same time.

**What was changed:**
* `rawstream.decode` works in blocks of 8 M words. The unwrapped clock is no longer
  built for every word, and there are no full-length index arrays. Output arrays are
  allocated once, after counting the record types. The per-photon word index
  (`ph_word`) was dropped; nothing used it.
* The TDC calibration is built from one joint (element × code) histogram instead of
  per-element copies.
* The pulse train is built in place in the memory of the sync time array, instead of
  5–6 full-length float64 arrays. The smoothing runs block-wise with the same
  arithmetic. The array is freed once photons and markers are located.
* Photon timing (fine time, delays, pulse, micro time) is computed in blocks of 4 M
  photons. The absolute time and phase are kept only when validating
  (`keep_validation_data`).
* Photons are sorted by pulse once, with a stable sort. Each output file then only
  inserts the few markers, instead of a full `lexsort` of all records per file. Records
  are encoded and written block-wise (`ptu_writer.write_ptu_t3`).
* The line-pass index is int32. It was not worth changing the word indices to int32,
  because they are no longer full-length.

**Result** (`/usr/bin/time -v`, peak resident memory, conversion without `--validate`):

| file | before | after |
|---|---|---|
| Matrix.raw (90 MB) | 1.21 GB, 2.2 s | 0.60 GB, 1.9 s |
| Sequence.raw (600 MB) | 7.19 GB, 27.8 s | 2.89 GB, 25.4 s |

**Acceptance — met:**
* PTU files byte-identical to the previous output, apart from `File_GUID` and
  `File_CreatingTime`. Checked for Matrix (both modes), Sequence (both modes,
  combined and split per sequence), the TDC tables and the conversion reports.
* All validation checks pass: Matrix 23/23, Sequence 25/25.
* Unit tests pass.

### B. Keep the sync records out of memory  *(done 2026-10-16)*

**Planned:** a compact sync model, i.e. one linear fit per block of about 64 sync
records instead of one smoothed time per record.

**Done instead (same goal, exact):** a linear fit per block would have changed the
micro times slightly. But the smoothing is local, so the existing smoothed sync times
can be **recomputed on demand** from the `.raw` file, block by block, with the same
arithmetic. The output stays bit-identical, and no per-record sync array is ever held
in memory.
* `rawstream.decode(..., keep_sync=False)` reads the file in blocks of 2 M words and
  never holds the whole `.raw`. It stores the unwrap state per block and the sync code
  histogram. `SyncReader` re-reads sync records for any range on demand.
* `timing.StreamedPulseTrain` (built by `streamed_pulse_train`) does three things:
  * one scan pass for the interval check, the missing records and the jitter;
  * the smoothed time at each block start, from values captured during that scan;
  * full 4 M-record blocks computed when photons need them (3 blocks cached).
  The moving mean is computed in sub-chunks with the running sum carried over, so
  the additions happen in the same order as one `np.cumsum`. Results are
  bit-identical to the in-memory `pulse_train_from_times`, which is still used for
  tiny files, by the tests and by `tools/make_figures.py`.
* Scan markers are located in the same sweep as the photons, so each pulse-train block
  is computed about once. Without this the run took 64 s instead of 31 s.
* The TDC calibration for the sync channel comes from the histogram. The sync-code
  histogram had been the 2.9 GB peak after A: a full-length int64 temporary.
* Smaller processing blocks: 2 M words, 1 M photons. Results do not depend on them;
  tests force tiny blocks.

**Result** (peak resident memory, conversion without `--validate`):

| file | before A | after A | after B |
|---|---|---|---|
| Matrix.raw (90 MB) | 1.21 GB, 2.2 s | 0.60 GB, 1.9 s | **0.44 GB**, 2.9 s |
| Sequence.raw (600 MB) | 7.19 GB, 27.8 s | 2.89 GB, 25.4 s | **1.49 GB**, 31.2 s |

The remaining peak is mostly per-photon data: the decoded photon records (11 bytes per
photon) plus the outputs (11 bytes per photon) during the photon sweep. On top of that
come a few hundred MB of pulse-train blocks and temporaries.

**Acceptance — met:** the PTU records, TDC tables, conversion reports and figures are
byte-identical to the previous version. Checked for Matrix and Sequence, both modes,
combined and split per sequence. Validation: 23/23 and 25/25. Unit tests: 12/12,
including streamed vs in-memory pulse train with gaps and tiny blocks, and chunked
vs one-shot smoothing.

## Open

### D. Memory of `--validate`

With validation, Sequence.raw still peaks at 7.3 GB, against 1.5 GB for the conversion alone. Validation decodes
the whole PTU with ptufile (twice: resolved and pooled), sorts it for the event
round trip and builds images at full size. Options: validate in blocks of records;
compare sorted per-pulse digests instead of full lexsorts; free the decoded records
between files.

### Discussed, not scheduled

* **C. Fully streaming conversion:** constant memory regardless of file size. After B,
  the `.raw` and the sync records are already streamed. What remains is per photon:
  * read photons from the file during the sweep instead of decoding them up front
    (saves 11 bytes per photon);
  * write the PTU without the global pulse sort, using a bounded reorder buffer.
  Only needed for acquisitions with many hundreds of millions of photons.
