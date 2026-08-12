# Benchmark and accuracy results

This directory keeps **raw measurements and their provenance**, not just screenshots or headline numbers.

## Datasets

### `pr4-fftw-baseline/`

The first controlled external production-library comparison. It compares fftlab's reusable complex/real plans with FFTW 3.3.10 under matched setup/execution semantics and two FFTW planner policies (`ESTIMATE`, `MEASURE`).

It contains:

- `metadata.json` — exact source commit/blobs, binary hash, FFTW runtime, compiler/build, protocol, and raw-file hashes;
- `raw/timings-session*.csv.gz` — every observation from three independently randomized sessions;
- `ANALYSIS.md` — execution speedups, bootstrap confidence intervals, setup economics, effect sizes, and amortized winners;
- `README.md` — the benchmark contract and interpretation scope.

The corpus contains **3,456 observations** across 6 sizes. Planning is separate from persistent-plan execution, FFTW inverse normalization is included in its execution timing, and specialized real FFTs are compared with specialized real FFTs.

```bash
python3 tools/analyze_vendor.py results/pr4-fftw-baseline/raw --seed 20260812
```

### `pr3-planning-real-baseline/`

The first formal dataset separating reusable FFT setup from steady-state execution and studying a specialized real-input representation. It contains **2,790 timing observations** across six power-of-two sizes and five modes.

```bash
python3 tools/analyze_plan.py results/pr3-planning-real-baseline/raw --seed 20260812
```

### `pr2-research-baseline/`

The first formal expanded-algorithm dataset: raw randomized timing samples, structural complexity models, multi-family forward/backward numerical error, bootstrap uncertainty, and effect-size analysis.

### `baseline-linux-amd-epyc-gcc14.csv`

The v1 development baseline, retained for historical continuity.

## Interpretation rules

1. Do not call an implementation universally “fastest” from these files.
2. Distinguish structural operation counts from measured latency.
3. Distinguish one-time planning/setup from reusable execution cost.
4. Do not substitute complex-input FFT measurements for real-input workloads.
5. Production-library comparisons must normalize transform convention, inverse scaling, planner policy, allocation/workspace semantics, and thread count.
6. Preserve every raw sample and exact source/runtime provenance.
7. Treat uncertainty intervals spanning parity as statistically unresolved rather than forcing a winner.
8. Update `SUMMARY.md` only with claims supported by checked-in data.

The current formal baselines are from a containerized/virtualized environment. They support implementation-specific and methodology claims, not universal hardware rankings.

Protocols: [`../docs/EXPERIMENTS.md`](../docs/EXPERIMENTS.md), [`../docs/PLANNING_REAL.md`](../docs/PLANNING_REAL.md), and [`../docs/VENDOR_BENCHMARKS.md`](../docs/VENDOR_BENCHMARKS.md).
