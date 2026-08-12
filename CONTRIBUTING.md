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

External production-library backends must additionally:

- remain optional to the from-first-principles core unless there is a documented reason otherwise;
- record exact runtime/library identity and version;
- pass deterministic numerical cross-checks before benchmark results are accepted;
- record planner flags/policy and cache/wisdom state;
- state how normalization differences are made equivalent;
- retain exact source, binary/build, and raw-data provenance;
- document alignment, allocator, threading, and build differences that remain construct-validity threats.

See [`docs/VENDOR_BENCHMARKS.md`](docs/VENDOR_BENCHMARKS.md).

## SIMD and architecture-specific kernels

Architecture-specific optimizations must not silently make the portable library require the optional ISA.

A SIMD/codelet contribution must:

- retain a scalar implementation or another documented portable fallback;
- gate every explicit optional ISA using a correct runtime capability check before execution;
- reject unsupported explicit ISA requests rather than executing undefined instructions;
- keep explicit kernel choices available even when an `Auto` policy exists, so dispatch can be audited;
- test each supported explicit path against an independent/reference implementation and verify inverse round trips;
- document compiler target requirements and the exact ISA/FMA feature set;
- report persistent-state and setup-cost changes introduced to improve vector access;
- separate algorithm changes from machine-code/layout changes wherever possible;
- preserve all explicit-width benchmark samples rather than reporting only the automatically selected winner;
- report auto-tuning cost and the selected policy per independent session;
- treat a wider vector width as a hypothesis, not an assumed win.

FMA kernels may differ from scalar code by floating-point rounding because fusion removes an intermediate rounding step. Correctness should therefore use a justified numerical tolerance rather than requiring bitwise equality unless bitwise reproducibility is itself the research objective.

An adaptive policy that occasionally disagrees with the best pooled explicit result is evidence about the tuner. Do not hard-code a benchmark-derived crossover merely to make the auto policy appear perfect. A static crossover may be introduced only with separately justified external-validity evidence and an explicit fallback/re-tuning strategy.

See [`docs/SIMD_KERNELS.md`](docs/SIMD_KERNELS.md).

## Benchmark changes

Do not optimize a benchmark to produce a preferred ranking. Preserve raw measurements, report the exact environment, and distinguish setup/planning from execution. Any change that alters timing semantics must update `docs/EXPERIMENTS.md`, the applicable specialized methodology, and result metadata.

For real-input studies, do not silently benchmark an N-point complex transform with zero imaginary input when the competing API exploits Hermitian symmetry.

For ISA studies, compare explicit implementations in the same randomized sessions whenever practical. Record CPU model, virtualization status, compiler, flags, capability set, and relevant library/runtime build. A CPUID feature bit is not performance evidence.

## Research claims

Claims should be specific enough to falsify. Prefer:

> On CPU X with compiler Y and flags Z, algorithm A had a median latency 1.24× lower than B at N=4096 across 5 sessions; bootstrap 95% speedup CI [1.19×, 1.29×].

Avoid:

> A is the fastest FFT.

For planning claims, state both setup cost and timed execution semantics. For vendor comparisons, report steady-state execution and setup/amortization when planner costs differ materially. For SIMD claims, state the explicit ISA, size range, compiler/environment, and uncertainty; do not generalize one measured AVX2/AVX-512 crossover to all processors.

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

Plan changes:

```bash
./build/fft-plan --self-test
```

FFTW adapter changes:

```bash
./build/fft-vendor --info
./build/fft-vendor --self-test
```

SIMD/codelet changes:

```bash
./build/fft-kernel --info
./build/fft-kernel --self-test
```

When multiple compilers are supported, architecture-specific code should be validated under each available compiler and under sanitizers where the sanitizer/toolchain supports the selected ISA paths.
