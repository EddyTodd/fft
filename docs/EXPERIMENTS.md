# Empirical research protocol

## Objective

Measure FFT implementations in a way that permits reproducible claims about **latency, variability, numerical accuracy, and crossover behavior**, while keeping theoretical arithmetic complexity separate from observed machine performance.

## Pre-registered hypotheses for the current implementation

H1. The direct DFT will show the expected quadratic scaling and rapidly lose to FFT families as N grows.

H2. Lower structural arithmetic count will not imply lower latency. In particular, this pedagogical recursive split-radix implementation may lose to iterative radix-2 because of allocation, recursion, trigonometric evaluation, and locality.

H3. Stockham’s regular autosort memory access will trade additional memory traffic/workspace for removal of explicit bit reversal; the crossover is architecture- and size-dependent.

H4. For prime sizes, Rader and Bluestein will exhibit measurable crossovers governed partly by convolution length and permutation/setup overhead.

H5. Numerical error will vary by algorithm and signal family even when all implementations pass conventional round-trip tests.

These are hypotheses, not README conclusions. Results should be updated only from checked-in raw data.

## Variables

### Independent variables

- algorithm;
- transform size N and factorization class;
- input signal family;
- compiler and version;
- optimization flags (`FFT_NATIVE`, LTO, build type);
- architecture/microarchitecture;
- session/order;
- optionally power mode, affinity, and thermal state.

### Response variables

- ns/transform for every raw timing sample;
- median latency;
- median absolute deviation (MAD);
- p05/p95;
- nonparametric 95% bootstrap CI for median;
- pairwise speedup with bootstrap CI;
- common-language probability of superiority (`P(faster)`);
- normalized forward/backward L1, L2, and Linf error;
- forward→inverse max absolute error;
- structural complex-operation and workspace models.

## Formal timing procedure

For publication-quality runs, use `tools/run_experiment.py` rather than only `--benchmark-suite`.

1. Build a release binary and record the exact commit/compiler/options.
2. Select sizes **before** looking at results. Include powers of two, smooth composites, awkward composites, and primes.
3. Run at least three independent sessions when making comparative claims.
4. Randomize `(N, algorithm)` execution order within each session to reduce monotonic thermal/time drift.
5. Perform untimed warmups.
6. Calibrate repetitions so each timing sample has sufficient duration.
7. Retain every raw sample; never keep only the best run.
8. Analyze medians and robust dispersion. Report confidence intervals and effect sizes, not only point estimates.
9. Treat overlapping intervals as descriptive uncertainty, not a mechanical hypothesis-test rule.
10. Re-run on materially different hardware before making architecture-general claims.

Example:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DFFT_NATIVE=ON
cmake --build build -j
python3 tools/run_experiment.py --binary build/fft --sessions 5 --samples 51
python3 tools/analyze.py results/run-*/timings.csv > analysis.md
```

## Timing semantics

The current transform functions return vectors and several pedagogical algorithms allocate scratch storage on each call. Therefore current timing is **end-to-end call latency including those allocations**. It is not equivalent to the execution-only timing of a production FFT plan with preallocated scratch/twiddle tables.

This is deliberate and must be reported. A future planner milestone should measure at least:

- plan/setup time;
- execution-only time with reused plan;
- amortized total cost for 1, 10, 100, and 10,000 transforms;
- in-place and out-of-place variants separately.

## Accuracy protocol

`--accuracy-suite` uses five deterministic input families:

- pseudorandom complex values;
- multi-tone signals plus noise;
- impulse;
- alternating/cancellation-heavy values;
- very high dynamic range values.

For each supported `(algorithm, N, signal)` tuple, compute:

1. long-double direct DFT reference;
2. normalized forward L1/L2/Linf error;
3. long-double inverse of the algorithm output to estimate backward error;
4. algorithmic forward→inverse round-trip max error.

The same generated signal is used for all algorithms at a given N and signal family.

## Threats to validity

### Internal

- OS scheduling and virtualization noise;
- frequency scaling/turbo/thermal throttling;
- benchmark process migration;
- calibration differences between very fast and very slow algorithms;
- allocation and page-fault effects;
- bootstrap samples are not proof of independence.

### Construct

- `5 N log2 N / time` is a radix-2-equivalent throughput scale, **not** measured FLOPs;
- complex-operation models do not equal instruction counts;
- long double is not arbitrary precision;
- end-to-end pedagogical call latency is not a vendor-library execution-plan benchmark.

### External

- one CPU/compiler cannot justify universal algorithm rankings;
- GPU, SIMD-specialized, multithreaded, real-input, multidimensional, and batched transforms can have qualitatively different winners.

## Publication rule

A claim belongs in `results/SUMMARY.md` only when all of the following are available:

- raw data;
- environment metadata;
- exact commit;
- command/protocol;
- sample/session counts;
- uncertainty/effect-size analysis;
- a scope statement saying which hardware/software configuration the claim applies to.
