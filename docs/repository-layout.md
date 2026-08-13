# Repository layout

`fftlab` follows the shared subject-repository contract defined by `EddyTodd/bench`.

## Permanent library

- `include/fftlab/` — installed public algorithm, plan, oracle, and kernel API.
- `src/kernel.cpp` — permanent non-header implementation.
- `tests/` — deterministic correctness tests.
- `docs/` — API, theory, planning, correctness, references, scope, and migration documentation.
- `cmake/` — package configuration.

Most FFT mechanisms are intentionally implemented in templated public headers, so a small permanent `src/` directory is expected.

## Transitional research

- `research/apps/` — historical benchmark/research executables.
- `research/legacy-src/` — development-era implementation/support sources no longer compiled by v1.
- `tools/` — historical Python campaign/analysis programs awaiting `bench` migration.
- `results/` — retained historical evidence.

The research relocation does not change the content of historical C++ sources. It prevents reviewers from mistaking uncompiled development artifacts for the active library implementation.
