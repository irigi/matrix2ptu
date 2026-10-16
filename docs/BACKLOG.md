# Backlog

Planned work that is agreed but not scheduled. See `REPORT.md` §9 for open questions.

## A. Reduce memory use: smaller and fewer arrays  *(agreed, not started)*

**Problem.** Peak memory is about 12× the `.raw` size (Sequence.raw, 600 MB: 7.2 GB for
the conversion, 8.5 GB with `--validate`). Most of it is not photons: sync records are
91 % (Matrix) and 77 % (Sequence) of all words, and several full-length 8-byte arrays
are kept at the same time.

**Changes:**
* Free intermediates as soon as they are used. For example, photon time and phase
  are not needed once `dtime` and the pulse index are computed.
* Use 32-bit integers for word indices and the line-pass index. This is valid for files
  below 2³¹ words, about 8 GB; check that limit and fall back to 64-bit above it.
* Keep the original record index (`raw_word`) only when validating.
* Replace the full `lexsort` of all output records in `_assemble` with a cheaper sort or
  merge. Photons are already almost in time order and markers are few.

**Expected:** peak about 5–6× the `.raw` size (≈ 3–4 GB for Sequence.raw).

**Acceptance:** PTU files byte-identical to the current output for Matrix.raw and
Sequence.raw. All validation checks pass. Measure peak memory before and after
(`/usr/bin/time -v`).

## Discussed, not scheduled

* **B. Compact sync model:** store one linear fit per block of about 64 sync records
  instead of one smoothed time per record. Must reproduce the current micro times.
* **C. Chunked (streaming) processing:** constant memory regardless of file size. Only
  needed for multi-GB acquisitions.
