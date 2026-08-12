# Contributing research and algorithms

Contributions are welcome when they preserve the repository’s two simultaneous goals: readable algorithm implementations and defensible empirical research.

## Algorithm changes

A new FFT implementation should include:

- a precise transform convention and supported-size domain;
- forward and inverse support where meaningful;
- deterministic correctness tests against an independent reference;
- at least one nontrivial numerical-accuracy case;
- a documented asymptotic complexity and workspace model;
- a primary literature reference when the algorithm is established work;
- inclusion in the benchmark suite only when comparisons are semantically fair.

## Benchmark changes

Do not optimize a benchmark to produce a preferred ranking. Preserve raw measurements, report the exact environment, and distinguish setup/planning from execution. Any change that alters timing semantics must update `docs/EXPERIMENTS.md` and result metadata.

## Research claims

Claims should be specific enough to falsify. Prefer:

> On CPU X with compiler Y and flags Z, algorithm A had a median latency 1.24x lower than B at N=4096 across 5 sessions; bootstrap 95% speedup CI [1.19x, 1.29x].

Avoid:

> A is the fastest FFT.

## Validation

Before opening a PR:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure

cmake -S . -B build-asan -DCMAKE_BUILD_TYPE=Debug -DFFT_SANITIZERS=ON
cmake --build build-asan -j
ctest --test-dir build-asan --output-on-failure
```
