# Research results summary

## v2 formal baseline

Environment: AMD EPYC 9V74 model exposed to a Linux container, GCC 14.2 release build, 3 independently randomized sessions, 31 timing samples per `(algorithm, N)` per session. Timings include per-call allocations. See `pr2-research-baseline/metadata.json` for full protocol metadata and exact source provenance.

### Timing

The fastest two paths at several sizes are often the `auto` wrapper and the exact implementation it dispatches to. Their tiny differences should not be interpreted as algorithmic wins when bootstrap speedup intervals span parity.

| N | Lowest measured median | Median | Comparator | Comparator / winner | Interpretation |
|---:|---|---:|---|---:|---|
| 64 | `auto` | 684.8 ns | radix2-iterative | 1.010× | unresolved; 95% speedup CI includes 1× |
| 127 | `rader` | 13.08 µs | auto | 1.005× | unresolved; both use the Rader path |
| 256 | `auto` | 3.220 µs | radix2-iterative | 1.002× | unresolved; same radix-2 path |
| 509 | `auto` | 59.29 µs | rader | 1.005× | unresolved; both use the Rader path |
| 1009 | `rader` | 130.01 µs | auto | 1.031× | small measured wrapper/path difference |
| 1024 | `radix2-iterative` | 15.66 µs | auto | 1.008× | unresolved; 95% speedup CI includes 1× |

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

This baseline is from a virtualized/containerized environment. It validates the research pipeline and supports implementation-specific hypotheses; it is **not** a universal hardware ranking. Stronger external-validity claims require repeated runs on named physical machines, multiple compilers, and controlled power/thermal/affinity conditions.

## v1 historical baseline

The original `baseline-linux-amd-epyc-gcc14.csv` remains available to show the development history of the benchmark harness. It used one 31-sample batch per algorithm/size and a smaller algorithm set. New research claims should prefer the v2 formal dataset because it preserves per-sample data, randomized multi-session order, confidence intervals, effect sizes, accuracy families, and exact source provenance.
