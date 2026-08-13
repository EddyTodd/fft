# Vendor benchmark status

Vendor-library timing is not part of fftlab's stable v1 library surface.

The development repository contains an FFTW adapter and historical controlled timing results. They were useful for understanding planning and implementation gaps, but vendor benchmarking, planner-effort comparisons, statistics, campaign orchestration, and cross-machine provenance belong in `EddyTodd/bench`.

The core `fftlab::fftlab` target has **no required FFTW dependency**. No vendor backend is required to build, install, or use v1 algorithms/plans.

When vendor research is migrated, preserve FFT-specific semantic normalization:

- forward sign and inverse normalization;
- real versus complex representation;
- in-place/out-of-place behavior;
- planning/setup versus reused execution;
- precision;
- thread count;
- allocation/alignment/workspace semantics.

See [`LEGACY_RESEARCH.md`](LEGACY_RESEARCH.md) for the files awaiting migration.
