# Contributing research and algorithms

Contributions are welcome when they preserve the repository’s two simultaneous goals: readable algorithm implementations and defensible empirical research.

## Algorithm changes

A new FFT implementation should include:

- a precise transform convention and supported-size domain;
- forward and inverse support where meaningful;
- deterministic correctness tests against an independent reference;
- at least one nontrivial numerical-accuracy case;
- documented asymptotic complexity and workspace/persistent-state models;
- a primary literature reference when the algorithm is established work;
- inclusion in benchmark suites only when comparisons are semantically fair.

## Plan and backend changes

Reusable plans and external-library backends must make lifecycle semantics explicit. Document and test:

- work performed during construction/setup;
- persistent state retained by the plan;
- caller/workspace buffers required by execution;
- whether execution allocates, performs trigonometric setup, or copies inputs;
- in-place/out-of-place behavior;
- complex versus real representation and output layout;
- normalization convention;
- thread count, precision, alignment, and ISA requirements where applicable.

Do not compare setup-inclusive latency from one implementation with reused-plan execution from another and label the result a general speed ranking.

## Arbitrary-length and prime plans

Changes to Bluestein, Rader, or a general arbitrary-length planner must make the **convolution reduction itself auditable**. Document:

- the convolution length chosen by each reduction;
- whether the convolution is direct cyclic or zero-padded linear-and-folded;
- precomputed chirp/kernel/permutation state;
- persistent memory and caller scratch requirements;
- setup and execution separately;
- explicit algorithm controls even when an automatic policy exists.

A prime dispatcher must not be justified only by one observed size threshold. Prefer structural rules that can be falsified across environments. If benchmark evidence changes a prior dispatcher, preserve the earlier evidence and explain why the semantics or implementation layer changed the conclusion.

When comparing Rader and Bluestein, setup-inclusive legacy calls may be used to quantify API/planning benefit, but **reusable algorithm ranking must use planned-vs-planned execution**.

See [`docs/ARBITRARY_PLANS.md`](docs/ARBITRARY_PLANS.md).

## External production-library backends

External backends must additionally:

- remain optional to the from-first-principles core unless documented otherwise;
- record exact runtime/library identity and version;
- pass deterministic numerical cross-checks before benchmark results are accepted;
- record planner flags/policy and cache/wisdom state;
- state how normalization differences are made equivalent;
- retain exact source, binary/build, and raw-data provenance;
- document alignment, allocator, threading, and build differences that remain validity threats.

See [`docs/VENDOR_BENCHMARKS.md`](docs/VENDOR_BENCHMARKS.md).

## SIMD and architecture-specific kernels

Architecture-specific optimizations must not silently make the portable library require the optional ISA.

A SIMD/codelet contribution must:

- retain a scalar implementation or another documented portable fallback;
- runtime-gate every optional ISA;
- reject unsupported explicit ISA requests;
- keep explicit kernel choices available even when `Auto` exists;
- test each supported explicit path against a reference and verify inverse round trips;
- document compiler target requirements and exact ISA/FMA feature sets;
- report persistent-state and setup-cost changes introduced for vector access;
- separate algorithm changes from machine-code/layout changes wherever possible;
- preserve all explicit-width benchmark samples rather than reporting only an automatic winner;
- report auto-tuning cost and selected policy per independent session;
- treat wider vectors as a hypothesis, not an assumed win.

FMA kernels may differ from scalar code by floating-point rounding because fusion removes an intermediate rounding step. Use justified numerical tolerances rather than requiring bitwise identity unless bitwise reproducibility is itself the research objective.

An adaptive policy that disagrees with the best pooled explicit result is evidence about the tuner. Do not hard-code a benchmark-derived crossover merely to make the auto policy appear perfect.

See [`docs/SIMD_KERNELS.md`](docs/SIMD_KERNELS.md).

## Benchmark changes

Do not optimize a benchmark to produce a preferred ranking. Preserve raw measurements, report the exact environment, and distinguish setup/planning from execution. Any change that alters timing semantics must update `docs/EXPERIMENTS.md`, the applicable specialized methodology, and result metadata.

For real-input studies, do not silently benchmark an N-point complex transform with zero imaginary input when the competing API exploits Hermitian symmetry.

For ISA studies, compare explicit implementations in the same randomized sessions whenever practical. Record CPU model, virtualization status, compiler, flags, capability set, and relevant library/runtime build. A CPUID feature bit is not performance evidence.

If tooling constraints require a reversible raw-data transport encoding, record the canonical raw-file hash and make the analyzer reconstruct the representation deterministically. Never replace raw observations with aggregates merely because binary upload is inconvenient.

## Research claims

Claims should be specific enough to falsify. Prefer:

> On CPU X with compiler Y and flags Z, algorithm A had a median latency 1.24× lower than B at N=4096 across 5 sessions; bootstrap 95% speedup CI [1.19×, 1.29×].

Avoid:

> A is the fastest FFT.

For planning claims, state both setup cost and timed execution semantics. For prime reductions, state the convolution lengths and whether Rader used a direct cyclic FFT. For vendor comparisons, report steady-state execution and setup/amortization when planner costs differ materially. For SIMD claims, state the explicit ISA, size range, compiler/environment, and uncertainty.

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

Specialized changes should also run their direct self-tests:

```bash
./build/fft-plan --self-test
./build/fft-vendor --info
./build/fft-vendor --self-test
./build/fft-kernel --info
./build/fft-kernel --self-test
./build/fft-arbitrary --info
./build/fft-arbitrary --self-test
```

When multiple compilers are supported, architecture- or planner-specific code should be validated under each available compiler and under sanitizers where supported.
