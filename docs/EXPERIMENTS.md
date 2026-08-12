# Empirical research protocol

## Objective

Measure FFT implementations in a way that permits reproducible claims about **latency, variability, numerical accuracy, planning cost, amortization, and crossover behavior**, while keeping theoretical arithmetic complexity separate from observed machine performance.

## Pre-registered hypotheses for the current implementation

H1. The direct DFT will show the expected quadratic scaling and rapidly lose to FFT families as N grows.

H2. Lower structural arithmetic count will not imply lower latency. In particular, this pedagogical recursive split-radix implementation may lose to iterative radix-2 because of allocation, recursion, trigonometric evaluation, and locality.

H3. Stockham’s regular autosort memory access will trade additional memory traffic/workspace for removal of explicit bit reversal; the crossover is architecture- and size-dependent.

H4. For prime sizes, Rader and Bluestein will exhibit measurable crossovers governed partly by convolution length and permutation/setup overhead.

H5. Numerical error will vary by algorithm and signal family even when all implementations pass conventional round-trip tests.

H6. For repeated power-of-two transforms, precomputing radix-2 permutation/twiddle state will materially reduce steady-state execution latency, and the one-time plan cost will amortize after a small finite number of transforms.

H7. For real-valued power-of-two input, an N/2-complex reduction with packed Hermitian output will materially outperform applying the already-planned N-point complex transform to data whose imaginary components are known to be zero.

These are hypotheses, not README conclusions. Results should be updated only from checked-in raw data.

## Variables

### Independent variables

- algorithm/decomposition;
- execution semantics: setup-inclusive call, reusable plan construction, or execution-only reused plan;
- data domain: complex or real input;
- transform size N and factorization class;
- input signal family;
- compiler and version;
- optimization flags (`FFT_NATIVE`, LTO, build type);
- architecture/microarchitecture;
- session/order;
- optionally power mode, affinity, and thermal state.

### Response variables

- ns/transform for every raw timing sample;
- plan/setup latency;
- setup break-even transform count;
- median latency;
- median absolute deviation (MAD);
- p05/p95;
- nonparametric 95% bootstrap CI for median;
- pairwise speedup with bootstrap CI;
- common-language probability of superiority (`P(faster)`), where applicable;
- normalized forward/backward L1, L2, and Linf error;
- forward→inverse max absolute error;
- structural complex-operation and workspace models;
- packed real-spectrum size and persistent-plan state where relevant.

## General algorithm timing procedure

For publication-quality comparisons among the original algorithm implementations, use `tools/run_experiment.py` rather than only `--benchmark-suite`.

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

## Timing semantics: pedagogical algorithm calls

The original transform functions return vectors and several pedagogical algorithms allocate scratch storage on each call. Their benchmark therefore measures **end-to-end call latency including implementation-specific allocation/setup performed by that function**.

This metric remains useful for comparing those APIs as implemented, but it must not be silently compared with an execution-only production FFT plan.

## Timing semantics: reusable plans

The v3 plan experiment resolves that ambiguity for reusable power-of-two radix-2 transforms. Its full contract and real-input derivation are in [`PLANNING_REAL.md`](PLANNING_REAL.md).

`tools/run_plan_experiment.py` / `fft-plan` distinguish five modes:

- complex plan construction;
- real plan construction;
- legacy radix-2 complex execution;
- planned radix-2 complex execution;
- planned real-input execution.

For the two planned execution modes:

- plans already exist before the timed interval;
- required output/scratch buffers already exist;
- execution performs no dynamic allocation;
- reusable trigonometric/twiddle setup is outside the timed interval;
- forward+inverse pairs are used so buffers can be reused without copying fresh input into every timed transform;
- pair elapsed time is divided by two;
- benchmark-mode order is randomized independently inside each sample;
- formal size order is randomized independently by session.

Setup is not discarded: it is measured independently and used to compute break-even transform counts.

Example:

```bash
python3 tools/run_plan_experiment.py \
  --binary build/fft-plan \
  --out results/run-plan \
  --sizes 64,256,1024,4096,16384,65536 \
  --sessions 5 \
  --samples 51 \
  --target-ms 5 \
  --source-commit "$(git rev-parse HEAD)"

python3 tools/analyze_plan.py results/run-plan/timings.csv
```

The analyzer accepts a single CSV, a gzip-compressed CSV, or a directory of CSV/CSV.GZ shards.

## Rules for future external-library comparisons

A production-library benchmark must normalize workload semantics before reporting rankings. Record separately whenever the API permits:

1. plan/setup construction time;
2. persistent plan memory;
3. caller-provided/workspace memory;
4. execution time using an already-created plan and already-allocated buffers;
5. destruction cost if the research question concerns one-shot transforms.

Also control and report:

- in-place vs out-of-place execution;
- real vs complex input;
- forward/inverse normalization convention;
- planning flags or effort levels;
- thread count;
- precision;
- alignment requirements;
- batch count;
- transform size/factorization;
- cold-plan latency vs repeated steady-state throughput.

A comparison that charges one implementation for planning on every transform while another reuses a plan is not an algorithm/library speed comparison unless that asymmetry exactly matches the intended application workload.

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

The real-input plan has an additional structural correctness check: every returned half-spectrum value is compared with the corresponding nonredundant bin of a full complex FFT of the same real sequence, followed by an inverse round trip.

## Threats to validity

### Internal

- OS scheduling and virtualization noise;
- frequency scaling/turbo/thermal throttling;
- benchmark process migration;
- calibration differences between very fast and very slow algorithms;
- allocation and page-fault effects;
- cache state and mode-order effects;
- bootstrap samples are not proof of independence.

### Construct

- `5 N log2 N / time` is a radix-2-equivalent throughput scale, **not** measured FLOPs;
- complex-operation models do not equal instruction counts;
- long double is not arbitrary precision;
- setup-inclusive pedagogical call latency is not equivalent to reused-plan execution latency;
- complex-input timing is not a faithful proxy for real-input workloads that can exploit Hermitian symmetry;
- break-even counts depend on both the specific plan constructor and the measured execution delta.

### External

- one CPU/compiler cannot justify universal algorithm rankings;
- GPU, SIMD-specialized, multithreaded, multidimensional, batched, and different-precision transforms can have qualitatively different winners;
- a virtualized/containerized development host is appropriate for validating the pipeline and local hypotheses, not for universal physical-hardware claims.

## Publication rule

A claim belongs in `results/SUMMARY.md` only when all of the following are available:

- raw data;
- environment metadata;
- exact source commit;
- command/protocol;
- sample/session counts;
- uncertainty/effect-size analysis;
- a scope statement saying which hardware/software configuration the claim applies to.

For reusable-plan or external-library claims, the publication record must additionally state whether planning, allocation, workspace preparation, data copying, and normalization are inside or outside the timed region.
