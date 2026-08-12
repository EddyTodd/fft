# Production-library benchmark contract

External FFT libraries are useful research controls only when the compared workloads are semantically equivalent. This document defines that contract and records the first backend: FFTW.

## FFTW adapter

`fft-vendor` loads FFTW dynamically at runtime. The core repository therefore has no build-time dependency on FFTW headers or libraries. If a compatible runtime is absent, the vendor self-test exits with code 77 and CTest records the test as skipped.

The adapter currently loads the double-precision serial API and benchmarks:

- complex in-place DFTs;
- real-to-complex / complex-to-real transforms with `N/2+1` nonredundant bins;
- `FFTW_ESTIMATE` planning;
- `FFTW_MEASURE` planning.

The runtime version and loaded library name are written into experiment metadata.

## Mathematical normalization

fftlab normalizes inverse transforms. FFTW's DFT convention leaves the inverse unnormalized. The benchmark therefore performs the required `1/N` multiplication **inside FFTW execution timing** after each inverse. This avoids crediting FFTW for omitting work required by the compared API contract.

The real comparison uses specialized real transforms on both sides. An N-point complex transform with zero imaginary input is not accepted as a proxy for a real FFT.

## Lifecycle normalization

Every strategy is decomposed into:

1. caller-buffer allocation;
2. plan/setup creation;
3. repeated execution with a persistent plan;
4. destruction/cleanup.

Caller-buffer allocation and destruction are outside setup/execution latency. Plan construction is measured separately. Execution uses preallocated buffers and already-created plans.

For FFTW, `MEASURE` setup is cold: wisdom is forgotten before every measured forward+inverse plan pair. This makes planner timing expensive by design but prevents earlier planning from contaminating later setup samples.

`ESTIMATE` and `MEASURE` execution use different persistent plan objects so their execution measurements correspond to the named planning policy.

## Timing procedure

For each transform size:

- generate deterministic complex and real input;
- create fftlab, FFTW ESTIMATE, and FFTW MEASURE plans;
- warm each execution path;
- calibrate a common iteration count;
- collect forward+inverse pair timings and divide by two;
- randomize all execution-mode order independently within every sample;
- collect cold setup observations separately;
- randomize setup-mode order;
- randomize size order independently across formal sessions;
- retain every observation.

The analyzer reports medians, bootstrap 95% speedup intervals, common-language effect sizes, setup break-even points, and amortized winners for representative repetition counts.

## What `ESTIMATE` and `MEASURE` mean

FFTW documents `FFTW_ESTIMATE` as a fast heuristic planner that avoids runtime measurements. `FFTW_MEASURE` actually measures candidate plans and can require substantially more planning time. The repository treats them as separate algorithms at the workload/lifecycle level rather than presenting a single ambiguous "FFTW" result.

The v4 formal run also calls `fftw_cleanup()` before unloading the runtime. This models FFTW's process-global planner lifecycle explicitly and keeps sanitizer validation meaningful.

## Alignment and implementation differences

FFTW buffers are allocated with `fftw_malloc`, while fftlab currently uses standard-library vector storage. Both meet their own backend requirements, but alignment and allocator behavior remain a construct-validity difference. A future controlled study should benchmark aligned fftlab buffers and explicitly inspect vectorization/code generation.

The current study is single-threaded. FFTW threading APIs are not enabled. It also does not compare wisdom persistence across processes, `PATIENT`/`EXHAUSTIVE`, arbitrary transform lengths, batches, multidimensional transforms, or SIMD-feature-matched custom FFTW builds.

## Backend admission requirements

A new production-library backend must:

- remain optional to the from-first-principles core;
- pass a deterministic numerical cross-check before timing;
- state transform convention and normalization;
- distinguish real from complex APIs;
- state in-place/out-of-place semantics;
- record planner policy and thread count;
- separate setup from reused execution;
- state allocation/alignment/workspace rules;
- retain raw observations and exact library/runtime provenance;
- avoid headline claims that exceed the measured hardware/software configuration.

Candidate next backends include Apple Accelerate/vDSP on macOS, Intel oneMKL on supported systems, and pocketfft as a compact modern CPU baseline.
