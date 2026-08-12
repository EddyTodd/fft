# Research results summary

## v3 planning and real-input baseline

The v3 experiment addresses a methodological flaw that matters for both this repository and future production-library comparisons: **one-time reusable setup must be separated from steady-state FFT execution**.

The checked-in `pr3-planning-real-baseline/` dataset contains 3 independently randomized sessions × 31 samples × 5 modes at each of 6 transform sizes, for **2,790 raw timing observations**. The benchmarked source is pinned to commit `372d59c592b9b7c8e30b7d168bf6c37449406509`; metadata also records the raw-data commit, file hashes, compiler/build configuration, and validation.

### Planning effect

`legacy-complex` and `planned-complex` implement the same radix-2 FFT decomposition. The difference is that the planned implementation precomputes bit reversal and twiddles once, then executes without allocation or trigonometric setup.

| N | Legacy complex | Planned complex | Speedup (bootstrap 95% CI) | Plan setup | Setup break-even |
|---:|---:|---:|---:|---:|---:|
| 64 | 661.3 ns | 415.2 ns | **1.593×** [1.578, 1.606] | 971 ns | 3.95 transforms |
| 256 | 3.27 µs | 2.12 µs | **1.542×** [1.522, 1.568] | 3.15 µs | 2.74 transforms |
| 1024 | 16.08 µs | 10.41 µs | **1.545×** [1.529, 1.582] | 11.14 µs | 1.96 transforms |
| 4096 | 77.65 µs | 51.65 µs | **1.504×** [1.482, 1.522] | 45.41 µs | 1.75 transforms |
| 16384 | 386.18 µs | 282.78 µs | **1.366×** [1.332, 1.406] | 213.38 µs | 2.06 transforms |
| 65536 | 1.775 ms | 1.370 ms | **1.295×** [1.273, 1.316] | 1.023 ms | 2.53 transforms |

Every speedup interval excludes parity. In this implementation and recorded environment, reusable planning reduces steady-state radix-2 latency by about **23–37%** relative to the legacy execution path, and plan construction repays itself after roughly **1.75–3.95 transforms**.

This changes how future library comparisons must be conducted. An FFTW/oneMKL/Accelerate execution benchmark should not be compared with a local implementation that reconstructs reusable twiddles on every call unless the research question explicitly includes setup in each operation.

### Real-input specialization

`RealRadix2Plan` exploits Hermitian symmetry. It packs an N-point real input into an N/2-point complex transform and returns only the N/2+1 nonredundant frequency bins.

| N | Planned complex | Planned real | Real-path speedup (bootstrap 95% CI) | Extra real setup break-even |
|---:|---:|---:|---:|---:|
| 64 | 415.2 ns | 350.1 ns | **1.186×** [1.177, 1.196] | 5.54 transforms |
| 256 | 2.12 µs | 1.61 µs | **1.320×** [1.298, 1.339] | 0.70 transforms |
| 1024 | 10.41 µs | 7.35 µs | **1.415×** [1.368, 1.438] | 0.22 transforms |
| 4096 | 51.65 µs | 32.72 µs | **1.579×** [1.558, 1.602] | 0.01 transforms |
| 16384 | 282.78 µs | 161.27 µs | **1.753×** [1.701, 1.804] | 0 transforms |
| 65536 | 1.370 ms | 785.55 µs | **1.744×** [1.708, 1.773] | 0 transforms |

The measured real-input advantage grows substantially with N in this baseline. For N≥256, any extra construction cost of the real plan is recovered within the first transform according to the recorded medians.

The broader conclusion is methodological rather than merely an optimization result: **complex-input FFT benchmarks are not faithful proxies for workloads whose data are known to be real**.

See `pr3-planning-real-baseline/ANALYSIS.md` and `docs/PLANNING_REAL.md` for raw-evidence interpretation and the mathematical/benchmark contract.

## v2 formal baseline

Environment: AMD EPYC 9V74 model exposed to a Linux container, GCC 14.2 release build, 3 independently randomized sessions, 31 timing samples per `(algorithm, N)` per session. Timings include per-call allocations. See `pr2-research-baseline/metadata.json` for full protocol metadata and exact source provenance.

### Timing

The fastest two paths at several sizes are often the `auto` wrapper and the exact implementation it dispatches to. Their tiny differences should not be interpreted as algorithmic wins when bootstrap speedup intervals span parity.

| N | Lowest measured median | Median | Runner-up / comparator | Interpretation |
|---:|---|---:|---|---|
| 64 | `auto` | 684.8 ns | radix2-iterative, 1.010× slower | unresolved; 95% speedup CI includes 1× |
| 127 | `rader` | 13.08 µs | auto, 1.005× slower | unresolved; both use the Rader path |
| 256 | `auto` | 3.220 µs | radix2-iterative, 1.002× slower | unresolved; same radix-2 path |
| 509 | `auto` | 59.29 µs | rader, 1.005× slower | unresolved; both use the Rader path |
| 1009 | `rader` | 130.01 µs | auto, 1.031× slower | small measured wrapper/path difference |
| 1024 | `radix2-iterative` | 15.66 µs | auto, 1.008× slower | unresolved; 95% speedup CI includes 1× |

The more informative cross-family comparisons are:

- N=509: Bluestein median 77.15 µs vs Rader 59.59 µs → **Rader ≈1.29× faster**.
- N=1009: Bluestein median 167.37 µs vs Rader 130.01 µs → **Rader ≈1.29× faster**.
- N=1024: direct DFT median 17.05 ms vs iterative radix-2 15.66 µs → **≈1089× speedup**.
- N=256: split-radix median 14.55 µs vs iterative radix-2 3.23 µs → split-radix is **≈4.51× slower** in this pedagogical implementation despite its lower structural multiplication model.

That last result is a central research point: arithmetic structure and machine performance are related but not interchangeable. Recursive allocation, data movement, twiddle evaluation, locality, and compiler behavior can dominate a theoretically attractive decomposition.

Use `pr2-research-baseline/ANALYSIS.md` for bootstrap intervals, MAD, p05/p95, complete rankings, and common-language effect sizes.

### Numerical accuracy

The formal accuracy subset uses a long-double direct DFT reference and five deterministic signal families.

At N=256, worst forward L2 error across the five inputs was:

| Algorithm | Worst forward L2 error |
|---|---:|
| split-radix | 2.97e-16 |
| radix-4 | 3.09e-16 |
| Stockham radix-2 | 3.52e-16 |
| mixed-radix | 7.01e-16 |
| iterative radix-2 | 2.56e-15 |
| Bluestein | 8.16e-15 |
| direct double DFT | 6.24e-14 |

The result does **not** imply a universal stability ranking. It demonstrates why floating-point accuracy needs an independent higher-precision reference and multiple signal families: directly evaluating the DFT formula in double precision is not automatically the most accurate computation.

At N=509, Rader and Bluestein are nearly tied in worst forward L2 error (~2.50e-14), while Rader is materially faster in this recorded environment. See `pr2-research-baseline/ACCURACY_ANALYSIS.md` for all rows.

## Scope of the evidence

The v2 and v3 formal baselines are from a virtualized/containerized environment. They validate the research pipeline and support implementation-specific hypotheses; they are **not** universal hardware rankings. Stronger external-validity claims require repeated runs on named physical machines, multiple compilers, and controlled power/thermal/affinity conditions.

## v1 historical baseline

The original `baseline-linux-amd-epyc-gcc14.csv` remains available to show the development history of the benchmark harness. It used one 31-sample batch per algorithm/size and a smaller algorithm set. New research claims should prefer the formal datasets because they preserve raw samples, randomized sessions, uncertainty analysis, scope statements, and exact provenance.
