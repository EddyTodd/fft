# Legacy empirical research notice

This file historically defined fftlab's benchmark-campaign protocol. Version 1.0 changes the permanent repository boundary: fftlab is an algorithm/planner library, while generic campaign orchestration, statistics, result collection, vendor timing, and provenance machinery are intended to move to `EddyTodd/bench`.

The historical experiment scripts and raw result corpora remain in this repository temporarily so evidence is not lost during the transition. They are **not part of the installed CMake package** and are not linked into `fftlab::fftlab`.

For the migration inventory see [`LEGACY_RESEARCH.md`](LEGACY_RESEARCH.md).

FFT-specific correctness methodology remains permanent:

- deterministic long-double direct DFT oracle;
- f32/f64 forward checks;
- normalized inverse round trips;
- transform identities;
- real Hermitian symmetry/layout validation;
- explicit planner/ISA mechanism cross-checks.

Generic benchmark methodology should be maintained in `bench` after migration rather than expanding this file again.
