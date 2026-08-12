# Empirical research protocol

## Objective

Measure FFT implementations in a way that permits reproducible claims about **latency, variability, numerical accuracy, planning cost, amortization, machine-level optimization, and crossover behavior**, while keeping mathematical complexity separate from observed implementation performance.

## Research layers

The repository intentionally has several benchmark layers because they answer different questions:

1. **algorithm calls** — direct DFT, radix families, mixed-radix, Rader, Bluestein, and dispatch as implemented;
2. **reusable plans** — construction separated from steady-state complex/real execution;
3. **production-library backends** — lifecycle-normalized external comparisons;
4. **machine-level codelets** — one mathematical decomposition held fixed while state layout and ISA change.

Do not combine results from these layers without stating the different timing semantics.

Specialized contracts:

- [`PLANNING_REAL.md`](PLANNING_REAL.md) — reusable plans and real-input reduction;
- [`VENDOR_BENCHMARKS.md`](VENDOR_BENCHMARKS.md) — external libraries;
- [`SIMD_KERNELS.md`](SIMD_KERNELS.md) — scalar/AVX2/AVX-512 codelets and plan-time ISA selection.

## Standing hypotheses

H1. Direct DFT will show quadratic scaling and rapidly lose to FFT families as N grows.

H2. Lower structural arithmetic count does not guarantee lower latency because allocation, recursion, locality, data movement, and code generation matter.

H3. Stockham's regular access trades workspace/traffic against explicit permutation; its crossover is architecture dependent.

H4. Rader and Bluestein exhibit prime-size crossovers governed by convolution length and setup/permutation effects.

H5. Numerical error varies by algorithm and signal family even when round-trip tests pass.

H6. Reusable radix-2 permutation/twiddle state materially reduces repeated execution latency and can amortize quickly.

H7. Specialized real FFTs materially outperform treating known-real data as generic complex data in relevant size ranges.

H8. SIMD codelets can materially close the gap to optimized production libraries without changing the mathematical FFT family.

H9. Wider SIMD is **not** assumed faster: AVX2/AVX-512 preference may cross over with N, microarchitecture, compiler, thermal/frequency state, and implementation layout.

H10. A plan-time empirical selector can adapt to hardware, but its own timing noise and construction cost may make it disagree with the best pooled explicit implementation. Such disagreements are results, not benchmark failures to erase.

Hypotheses become conclusions only from checked-in raw evidence.

## Variables

### Independent variables

- algorithm/decomposition;
- implementation/codelet and explicit ISA;
- planning/dispatch policy;
- setup-inclusive call versus reusable execution;
- complex versus real data;
- N and factorization class;
- signal family;
- compiler/version and optimization flags;
- architecture/microarchitecture and runtime capability set;
- vendor runtime/build and planner policy;
- session and randomized order;
- optionally affinity, power mode, thermal state, alignment, and thread count.

### Response variables

- every raw ns/transform observation;
- setup/plan/tuning latency;
- setup break-even count;
- median, MAD, p05/p95;
- bootstrap confidence intervals for medians and speedup ratios;
- effect-size/probability-of-superiority metrics where useful;
- forward/backward normalized L1/L2/L∞ error;
- round-trip error;
- structural operation/workspace models;
- persistent plan/codelet storage;
- selected policy/ISA for adaptive planners.

## General timing procedure

Publication-quality performance claims should satisfy the following unless a specialized protocol documents a justified exception:

1. Pin the exact source commit and record compiler/options.
2. Select the transform matrix before interpreting results.
3. Use at least three independent sessions for comparative claims.
4. Randomize competing mode order within samples and size/order across sessions where practical.
5. Perform untimed warmups.
6. Calibrate repeated operations so timer overhead is small relative to sample duration.
7. Retain **every** raw sample, not just best runs or aggregates.
8. Report medians and robust dispersion plus uncertainty/effect size.
9. Do not use interval overlap as a mechanical significance test; compute the relevant pairwise statistic directly.
10. Re-run on materially different physical hardware before architecture-general claims.

## Original algorithm calls

The historical algorithm APIs may allocate or generate reusable state on every call. Their benchmark is therefore end-to-end latency **for those functions as implemented**.

That metric is valid for the API but must not be silently compared to reused-plan execution from another implementation.

```bash
python3 tools/run_experiment.py --binary build/fft --sessions 5 --samples 51
python3 tools/analyze.py results/run-*/timings.csv
```

## Reusable plan timing

For planned execution:

- plans exist before execution timing;
- output/scratch buffers already exist;
- reusable trigonometric/twiddle setup is outside execution timing;
- forward+inverse pairs permit buffer reuse without timing a fresh input copy;
- elapsed pair time is divided by two;
- setup is measured separately and used in amortization.

```bash
python3 tools/run_plan_experiment.py \
  --binary build/fft-plan --out results/run-plan \
  --sizes 64,256,1024,4096,16384,65536 \
  --sessions 5 --samples 51 --target-ms 5 \
  --source-commit "$(git rev-parse HEAD)"
```

## Production-library timing

External comparisons must normalize, record, or explicitly identify remaining differences in:

- setup/planning policy and cache/wisdom state;
- persistent plan and workspace ownership;
- in-place/out-of-place semantics;
- complex/real representation;
- inverse normalization;
- precision;
- alignment/allocator requirements;
- threads and batch count;
- library identity/version/build;
- cold versus reused execution.

A faster persistent kernel can have worse total workload cost if planning is expensive. Both quantities should be reported when material.

## SIMD/codelet timing

ISA studies must distinguish **capability**, **explicit implementation performance**, and **adaptive selection**.

For the v5 radix-2 codelet study:

- scalar, AVX2/FMA, and AVX-512/FMA share the same swap-list and stage-contiguous twiddle plan structure;
- the mathematical radix-2 decomposition is held fixed;
- explicit-width modes remain available even when `Auto` is tested;
- the portable build does not globally require optional ISA support;
- runtime capability checks gate optional instruction paths;
- unsupported explicit ISA requests are rejected;
- FMA paths are checked numerically with a justified tolerance rather than bitwise equality;
- plan/tuning cost is measured separately;
- auto-selected ISA is retained per session;
- explicit kernels and auto are timed in the same randomized matrix as the v3 plan and FFTW when available.

Formal example:

```bash
python3 tools/run_kernel_experiment.py \
  --binary build/fft-kernel \
  --out results/run-kernel \
  --sizes 64,256,1024,4096,16384,65536 \
  --sessions 3 --samples 31 --setup-samples 1 \
  --target-ms 2 \
  --source-commit "$(git rev-parse HEAD)"

python3 tools/analyze_kernel.py results/run-kernel/raw \
  --bootstrap 5000 --seed 20260812
```

A wider vector implementation is called faster only from observed data with a stated uncertainty interval. Do not infer a universal ISA crossover from one machine.

An adaptive selector should be evaluated against the explicit paths it can choose. If it disagrees with the pooled explicit winner, preserve and report the disagreement. Do not tune the benchmark or hard-code the observed crossover solely to make the adaptive policy appear correct.

## Accuracy protocol

The general accuracy suite uses five deterministic inputs:

- pseudorandom complex values;
- multi-tone plus noise;
- impulse;
- alternating/cancellation-heavy values;
- high-dynamic-range values.

For each supported `(algorithm, N, signal)` tuple:

1. compute a long-double direct DFT reference;
2. report normalized forward L1/L2/L∞ error;
3. use a long-double inverse to estimate backward error;
4. report algorithmic forward→inverse round-trip max error.

Reusable real plans additionally verify each packed half-spectrum bin against the corresponding full complex FFT bin. SIMD kernels are cross-checked against the established radix-2 result and round-trip behavior; FMA-induced rounding differences are permitted within the declared tolerance.

## Threats to validity

### Internal

- OS scheduling and virtualization noise;
- frequency/turbo/thermal state;
- process migration and affinity;
- cache state and mode order;
- calibration differences;
- page faults/allocation outside or inside the intended timing region;
- adaptive tuner noise;
- bootstrap samples do not prove independence.

### Construct

- `5 N log2 N / time` is a radix-2-equivalent scale, not measured FLOPs;
- complex structural counts are not machine instruction counts;
- long double is not arbitrary precision;
- setup-inclusive calls and reused-plan execution are different workloads;
- complex timing is not a proxy for a specialized real FFT;
- ISA availability is not throughput evidence;
- FMA changes floating-point rounding semantics without changing the intended DFT;
- a different plan layout may change memory footprint as well as execution speed;
- an adaptive planner must be charged for tuning.

### External

- one virtualized CPU/compiler cannot justify universal algorithm or ISA rankings;
- physical machines may have different AVX-512 frequency behavior, cache hierarchy, memory system, compiler output, and scheduler state;
- Arm SIMD, GPUs, threads, batching, multidimensional transforms, and different precisions can have qualitatively different winners.

## Publication rule

A headline claim belongs in `results/SUMMARY.md` only when the repository contains:

- all raw data;
- environment/runtime metadata;
- exact source commit and relevant source blobs;
- build/compiler provenance;
- command/protocol and timing semantics;
- sample/session counts;
- uncertainty/effect-size analysis;
- an explicit scope statement.

Reusable/vendor claims must state planning, allocation/workspace, data-copy, normalization, and threading semantics. SIMD claims must additionally record capability requirements, explicit-width controls, adaptive selections if any, and persistent-state/setup tradeoffs.
