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

## Plan and backend changes

Reusable plans and external-library backends must make lifecycle semantics explicit. Document and test:

- what work occurs during plan/setup construction;
- what persistent state the plan retains;
- what caller/workspace buffers execution requires;
- whether execution allocates, performs trigonometric setup, or copies inputs;
- in-place/out-of-place behavior;
- complex vs real representation and output layout;
- normalization convention;
- thread count/precision/alignment requirements where applicable.

Do not compare setup-inclusive latency from one implementation with reused-plan execution from another and label the result a general speed ranking. If both workloads are important, measure both separately.

External production-library backends must additionally:

- remain optional to the from-first-principles core unless there is a documented reason otherwise;
- record the exact runtime/library identity and version;
- pass deterministic numerical cross-checks before benchmark results are accepted;
- record planner flags/policy and whether caches or wisdom were warm or cold;
- state how normalization differences are made equivalent;
- retain exact source, binary/build, and raw-data provenance;
- document alignment, allocator, and threading differences that remain construct-validity threats.

The current production-library contract is [`docs/VENDOR_BENCHMARKS.md`](docs/VENDOR_BENCHMARKS.md).

## Benchmark changes

Do not optimize a benchmark to produce a preferred ranking. Preserve raw measurements, report the exact environment, and distinguish setup/planning from execution. Any change that alters timing semantics must update `docs/EXPERIMENTS.md`, the applicable specialized methodology (`docs/PLANNING_REAL.md` and/or `docs/VENDOR_BENCHMARKS.md`), and result metadata.

For real-input studies, do not silently benchmark an N-point complex transform with zero imaginary input when the competing API exploits Hermitian symmetry. Treat complex and real transforms as distinct workloads.

## Research claims

Claims should be specific enough to falsify. Prefer:

> On CPU X with compiler Y and flags Z, algorithm A had a median latency 1.24x lower than B at N=4096 across 5 sessions; bootstrap 95% speedup CI [1.19x, 1.29x].

Avoid:

> A is the fastest FFT.

For planning claims, state both setup cost and the timed execution semantics. For example:

> With a reused N=4096 radix-2 plan, execution was 1.50x faster than the setup-inclusive legacy path; median plan construction amortized after 1.8 transforms in this environment.

For vendor comparisons, report both steady-state execution and setup/amortization when planner costs differ materially. A faster execution kernel is not automatically a cheaper workload for one-shot or short-run use.

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

Changes touching the plan layer should also run:

```bash
./build/fft-plan --self-test
```

Changes touching a vendor adapter should run that adapter's deterministic self-test on every available supported compiler/runtime configuration before publishing measurements. For FFTW:

```bash
./build/fft-vendor --info
./build/fft-vendor --self-test
```
