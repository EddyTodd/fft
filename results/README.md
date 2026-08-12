# Benchmark and accuracy results

This directory keeps **raw measurements and their provenance**, not just screenshots or headline numbers.

## Datasets

### `pr3-planning-real-baseline/`

The first formal dataset that separates reusable FFT setup from steady-state execution and studies a specialized real-input representation. It contains:

- `metadata.json` — exact benchmarked source commit, raw-data commit, environment, compiler/build, protocol, hashes, and validation;
- `raw/timings-session*.csv.gz` — every raw timing observation for three randomized sessions;
- `ANALYSIS.md` — bootstrap speedup intervals, setup medians, and plan-amortization estimates;
- `README.md` — benchmark-mode definitions and interpretation guidance.

The dataset contains **2,790 timing observations** across 6 power-of-two sizes and five modes: complex-plan construction, real-plan construction, legacy complex execution, planned complex execution, and planned real execution. Setup is measured separately rather than folded into repeated execution latency.

Reproduce its analysis directly from the compressed checked-in data:

```bash
python3 tools/analyze_plan.py results/pr3-planning-real-baseline/raw \
  --seed 20260812
```

### `pr2-research-baseline/`

The first formal research dataset for the expanded algorithm set. It contains:

- `metadata.json` — exact source commit, environment, and protocol;
- `raw/timings-N*-session*.csv` — every raw timing sample with session and randomized execution order;
- `complexity.csv` — structural operation/workspace models for measured algorithm/size pairs;
- `raw/accuracy-N*.csv` — forward/backward normed error across deterministic signal families;
- `ANALYSIS.md` — timing rankings, bootstrap confidence intervals, robust dispersion, and effect sizes;
- `ACCURACY_ANALYSIS.md` — accuracy ranking by worst forward L2 error across signal families.

### `baseline-linux-amd-epyc-gcc14.csv`

The original v1 development baseline. It predates the formal randomized multi-session protocol and is retained for historical continuity.

## Interpretation rules

1. Do not call an implementation universally “fastest” from these files.
2. Distinguish structural operation counts from measured latency.
3. Distinguish one-time planning/setup from reusable execution cost.
4. Do not substitute complex-input FFT measurements for real-input workloads when Hermitian symmetry can be exploited.
5. Preserve raw samples and metadata when adding a result.
6. Treat tiny differences with uncertainty intervals spanning parity as statistically unresolved.
7. Update `SUMMARY.md` only with claims supported by checked-in data.

The current datasets are from a containerized/virtualized environment, so absolute timings are **pipeline/research-development evidence**, not a substitute for results from named physical machines with controlled power, affinity, and thermal state.

The general protocol is in [`../docs/EXPERIMENTS.md`](../docs/EXPERIMENTS.md); setup/execution and real-input semantics are specified in [`../docs/PLANNING_REAL.md`](../docs/PLANNING_REAL.md).
