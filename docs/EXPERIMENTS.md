# Empirical research protocol

## Objective

Measure FFT implementations in a way that permits reproducible claims about **latency, variability, numerical accuracy, planning cost, amortization, machine-level optimization, reduction structure, and crossover behavior**, while keeping mathematical complexity separate from observed implementation performance.

## Research layers

The repository intentionally maintains distinct benchmark layers because they answer different questions:

1. **algorithm calls** — historical direct/radix/mixed/Rader/Bluestein APIs as implemented, including per-call setup/allocation where present;
2. **reusable power-of-two plans** — construction separated from steady-state complex/real execution;
3. **reusable arbitrary-length plans** — Bluestein/Rader convolution kernels and permutations made persistent;
4. **production-library backends** — lifecycle-normalized external comparisons;
5. **machine-level codelets** — mathematical decomposition held fixed while layout and ISA change.

Do not combine results from different layers without stating their timing semantics.

Specialized contracts:

- [`PLANNING_REAL.md`](PLANNING_REAL.md) — reusable radix-2 plans and real-input reduction;
- [`ARBITRARY_PLANS.md`](ARBITRARY_PLANS.md) — reusable Bluestein/Rader plans and prime dispatch;
- [`VENDOR_BENCHMARKS.md`](VENDOR_BENCHMARKS.md) — external libraries;
- [`SIMD_KERNELS.md`](SIMD_KERNELS.md) — scalar/AVX2/AVX-512 codelets and plan-time ISA selection.

## Standing hypotheses

H1. Direct DFT exhibits quadratic scaling and rapidly loses to FFT families as N grows.

H2. Lower structural arithmetic count does not guarantee lower latency because allocation, recursion, locality, data movement, and code generation matter.

H3. Stockham's regular access trades workspace/traffic against explicit permutation; its crossover is architecture dependent.

H4. Rader and Bluestein exhibit prime-size crossovers governed by convolution length and setup/permutation effects.

H5. Numerical error varies by algorithm and signal family even when round-trip tests pass.

H6. Reusable radix-2 permutation/twiddle state materially reduces repeated execution latency and can amortize quickly.

H7. Specialized real FFTs materially outperform treating known-real data as generic complex data in relevant size ranges.

H8. SIMD codelets can materially close the gap to optimized production libraries without changing the mathematical FFT family.

H9. Wider SIMD is not assumed faster: AVX2/AVX-512 preference may cross over with N, microarchitecture, compiler, frequency state, and implementation layout.

H10. A plan-time empirical ISA selector can adapt to hardware, but timing noise and construction cost may make it disagree with the best pooled explicit path. Such disagreements are results, not failures to erase.

H11. Persistent Bluestein/Rader state materially changes prime-transform execution cost relative to setup-inclusive historical APIs.

H12. Rader's strongest reusable advantage occurs when `p-1` permits a substantially shorter convolution than Bluestein, especially a direct power-of-two cyclic FFT.

H13. Equal convolution lengths do not imply a universal Rader/Bluestein winner because chirp products, permutations, folding, and memory traffic differ.

H14. A reusable prime dispatcher should prefer structural evidence over a machine-specific size threshold and remain auditable through explicit algorithm controls.

Hypotheses become conclusions only from checked-in evidence.

## Variables

### Independent variables

- algorithm/decomposition or convolution reduction;
- implementation/codelet and explicit ISA;
- planning/dispatch policy;
- setup-inclusive call versus reusable execution;
- complex versus real data;
- N and factorization class, including factorization of `p-1` for Rader;
- convolution length and direct-cyclic versus zero-padded mode;
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
- normalized forward/backward L1/L2/L∞ error;
- round-trip error;
- structural operation/workspace models;
- persistent plan/codelet storage;
- selected policy/ISA for adaptive planners;
- convolution length and reduction mode for arbitrary-length plans.

## General timing procedure

Publication-quality performance claims should satisfy the following unless a specialized protocol documents a justified exception:

1. Pin the exact source commit and record compiler/options.
2. Select the transform matrix before interpreting results.
3. Use at least three independent sessions for comparative claims.
4. Randomize competing mode order within samples and size order across sessions where practical.
5. Perform untimed warmups.
6. Calibrate repeated operations so timer overhead is small relative to sample duration.
7. Retain **every** raw sample, not just best runs or aggregates.
8. Report medians and robust dispersion plus direct pairwise uncertainty/effect size.
9. Do not use interval overlap as a mechanical significance test; compute the relevant pairwise statistic.
10. Re-run on materially different physical hardware before architecture-general claims.

If an execution environment/tool cannot upload a canonical binary raw artifact, a reversible transport encoding is acceptable only when the canonical raw hash is retained and analysis reconstructs the original observations deterministically.

## Historical algorithm calls

The original algorithm APIs may allocate or generate reusable state on every call. Their benchmark is end-to-end latency **for those functions as implemented**. This remains useful for API/history studies but is not a reusable-plan ranking.

```bash
python3 tools/run_experiment.py --binary build/fft --sessions 5 --samples 51
python3 tools/analyze.py results/run-*/timings.csv
```

## Reusable power-of-two plans

For planned execution:

- plans exist before execution timing;
- output/scratch buffers already exist;
- reusable trigonometric/twiddle setup is outside execution timing;
- forward+inverse pairs permit buffer reuse without timing fresh input copies;
- elapsed pair time is divided by two;
- setup is measured separately and used in amortization.

```bash
python3 tools/run_plan_experiment.py \
  --binary build/fft-plan --out results/run-plan \
  --sizes 64,256,1024,4096,16384,65536 \
  --sessions 5 --samples 51 --target-ms 5 \
  --source-commit "$(git rev-parse HEAD)"
```

## Reusable arbitrary-length plans

The v6 prime benchmark separates two questions:

1. how much planning improves each historical reduction;
2. which reduction is faster **after both are planned fairly**.

Execution modes:

- legacy Bluestein;
- planned Bluestein;
- legacy Rader;
- planned Rader;
- FFTW `ESTIMATE`;
- FFTW `MEASURE`.

The legacy modes remain setup-inclusive and are used only for planning-benefit estimates. Prime algorithm/dispatch conclusions use planned-vs-planned observations.

For planned Bluestein/Rader:

- chirps/kernel spectra/permutations/radix-2 convolution plans already exist;
- caller-owned output and scratch vectors already exist;
- no trigonometric setup or dynamic allocation occurs during planned execution;
- forward+inverse pairs are divided by two;
- inverse normalization remains inside the timed path;
- the harness records convolution M and whether Rader used direct cyclic convolution;
- all six execution modes share a common calibrated iteration count within each N;
- mode order is randomized per sample and size order per session.

Formal example:

```bash
python3 tools/run_arbitrary_experiment.py \
  --binary build/fft-arbitrary \
  --out results/run-arbitrary \
  --sizes 17,31,61,127,257,509,1009,4093 \
  --sessions 3 --samples 31 --setup-samples 1 --target-ms 2 \
  --source-commit "$(git rev-parse HEAD)"

python3 tools/analyze_arbitrary.py results/run-arbitrary/raw \
  --bootstrap 5000 --seed 20260812
```

A prime automatic rule must be evaluated against explicit planned Bluestein and Rader modes. A tiny one-host advantage when convolution lengths tie is insufficient by itself to justify a permanent static crossover. Structural rules should remain revisable after mixed-radix-convolution or physical-hardware studies.

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

ISA studies distinguish **capability**, **explicit implementation performance**, and **adaptive selection**.

For the v5 radix-2 study:

- scalar, AVX2/FMA, and AVX-512/FMA share matched plan structure;
- the radix-2 decomposition is held fixed;
- explicit-width modes remain available even when `Auto` is tested;
- portable builds do not globally require optional ISAs;
- runtime checks gate optional instructions;
- unsupported explicit ISA requests are rejected;
- FMA paths use justified numerical tolerance rather than bitwise identity;
- plan/tuning cost is measured separately;
- auto-selected ISA is retained per session;
- explicit kernels and auto run in the same randomized matrix as the v3 plan and FFTW.

A wider vector implementation is called faster only from observed data with stated uncertainty. An adaptive selector that disagrees with the pooled explicit winner remains part of the evidence.

## Accuracy protocol

The general accuracy suite uses deterministic pseudorandom, multi-tone/noise, impulse, alternating/cancellation-heavy, and high-dynamic-range inputs.

For each supported `(algorithm, N, signal)` tuple:

1. compute a long-double direct DFT reference;
2. report normalized forward L1/L2/L∞ error;
3. use a long-double inverse to estimate backward error;
4. report forward->inverse round-trip max error.

Reusable real plans additionally verify each packed half-spectrum bin against the corresponding full complex FFT bin. SIMD kernels are cross-checked against established radix-2 output. Arbitrary plans are directly compared with the DFT at bounded sizes, round-trip tested, and cross-checked against FFTW on selected primes.

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
- structural complex-operation counts are not machine instruction counts;
- long double is not arbitrary precision;
- setup-inclusive calls and reused-plan execution are different workloads;
- complex timing is not a proxy for specialized real FFTs;
- equal convolution lengths do not make Rader/Bluestein implementations identical;
- the current Rader fallback's radix-2 zero padding is part of the implementation under study, not a theoretical requirement;
- ISA availability is not throughput evidence;
- FMA changes rounding semantics without changing the intended DFT;
- different plan layouts change memory footprint as well as latency;
- adaptive planning must be charged for tuning.

### External

- one virtualized CPU/compiler cannot justify universal algorithm, prime-dispatch, or ISA rankings;
- physical machines may have different cache hierarchy, memory system, compiler output, vector-frequency behavior, and scheduler state;
- mixed-radix convolution, Arm SIMD, GPUs, threads, batching, multidimensional transforms, and other precisions can have qualitatively different winners.

## Publication rule

A headline claim belongs in `results/SUMMARY.md` only when the repository contains:

- all raw observations or a lossless/hash-verifiable transport of them;
- environment/runtime metadata;
- exact source commit and relevant source blobs;
- build/compiler provenance;
- command/protocol and timing semantics;
- sample/session counts;
- uncertainty/effect-size analysis;
- an explicit scope statement.

Reusable/vendor claims must state planning, allocation/workspace, data-copy, normalization, and threading semantics. Arbitrary-plan claims must state convolution lengths and Rader's cyclic/folded mode. SIMD claims must record capability requirements, explicit-width controls, adaptive selections if any, and persistent-state/setup tradeoffs.
