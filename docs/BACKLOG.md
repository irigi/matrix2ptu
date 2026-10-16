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

## Open

### D. Memory of `--validate`

With validation, Sequence.raw still peaks at 7.2 GB (was 8.5 GB). Validation decodes
the whole PTU with ptufile (twice: resolved and pooled), sorts it for the event
round trip and builds images at full size. Options: validate in blocks of records;
compare sorted per-pulse digests instead of full lexsorts; free the decoded records
between files.

### Discussed, not scheduled

* **B. Compact sync model:** store one linear fit per block of about 64 sync records
  instead of one smoothed time per record. After A, the sync array is the largest
  remaining item during conversion (8 bytes per sync record: 0.9 GB for Sequence).
  Must reproduce the current micro times.
* **C. Chunked (streaming) processing:** constant memory regardless of file size. Only
  needed for multi-GB acquisitions. Includes not reading the whole `.raw` into memory.
