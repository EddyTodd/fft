# v3 planning and real-input baseline

This result package measures the distinction introduced in v3 between **one-time plan construction** and **steady-state FFT execution**, and separately measures the benefit of exploiting real-input Hermitian symmetry.

## Evidence package

- [`metadata.json`](metadata.json) records environment, protocol, exact source/raw-data commits, validation, file hashes, and timing semantics.
- [`raw/`](raw/) contains all **2,790 raw timing observations** as three gzip-compressed CSV session files. Each file contains exactly 930 observations: 6 sizes × 5 benchmark modes × 31 samples.
- [`ANALYSIS.md`](ANALYSIS.md) is regenerated from the raw session files with `tools/analyze_plan.py`.

Reproduce the checked-in analysis after cloning:

```bash
python3 tools/analyze_plan.py results/pr3-planning-real-baseline/raw \
  --seed 20260812 \
  --output reproduced.md
```

## Benchmark modes

| Mode | Meaning |
|---|---|
| `complex-setup` | construct a fresh N-point `Radix2Plan` |
| `real-setup` | construct a fresh N-point `RealRadix2Plan` |
| `legacy-complex` | existing in-place radix-2 forward/inverse execution with repeated twiddle generation |
| `planned-complex` | reused N-point plan; precomputed permutation/twiddles; no allocation or trig setup in execution |
| `planned-real` | reused real-input plan; N/2 complex FFT plus half-spectrum recombination |

Execution samples time a forward+inverse pair and divide elapsed time by two. Mode order is randomized inside every sample; transform-size order is randomized independently for each session. Setup is never folded into execution latency.

## Main findings

Across every tested N, the planned complex path is materially faster than the mathematically equivalent legacy radix-2 path; all reported bootstrap 95% speedup intervals exclude 1×. Median speedup ranges from **1.295× to 1.593×**.

Complex-plan construction amortizes after approximately **1.75–3.95 transforms** over the tested sizes. This means repeated transforms reach the steady-state benefit very quickly in this implementation.

The real-input specialization is also faster than the already-planned complex path. The advantage increases with transform size in this baseline, from **1.186× at N=64** to **1.744× at N=65536**. For N≥256, the extra real-plan construction cost is repaid within the first transform according to the recorded medians.

These findings support two important methodological conclusions:

1. repeated-transform benchmarks should not mix reusable setup work into every execution unless that is explicitly the workload being studied;
2. complex FFT benchmarks are not substitutes for real-input FFT benchmarks when the application data is known to be real.

## Scope

The baseline was captured in a virtualized/containerized environment exposing an AMD EPYC 9V74 model, using GCC 14.2 with an `-O3 -DNDEBUG -std=c++23` build equivalent to `FFT_NATIVE=OFF` and no LTO. It is evidence about this implementation and the benchmark protocol, not a universal machine ranking.

The next production-library comparison should preserve the same setup/execution separation and real-vs-complex distinction. See [`docs/PLANNING_REAL.md`](../../docs/PLANNING_REAL.md).
