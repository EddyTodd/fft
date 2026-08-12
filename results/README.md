# Benchmark and accuracy results

This directory keeps **raw measurements and their provenance**, not just screenshots or headline numbers.

## Datasets

### `baseline-linux-amd-epyc-gcc14.csv`

The original v1 development baseline. It predates the formal randomized multi-session protocol and is retained for historical continuity.

### `pr2-research-baseline/`

The first formal research dataset for the expanded algorithm set. It contains:

- `metadata.json` — exact source commit, environment, and protocol;
- `raw/timings-N*-session*.csv` — every raw timing sample with session and randomized execution order;
- `complexity.csv` — structural operation/workspace models for measured algorithm/size pairs;
- `raw/accuracy-N*.csv` — forward/backward normed error across deterministic signal families;
- `ANALYSIS.md` — timing rankings, bootstrap confidence intervals, robust dispersion, and effect sizes;
- `ACCURACY_ANALYSIS.md` — accuracy ranking by worst forward L2 error across signal families.

The environment is containerized/virtualized, so absolute timings are **pipeline/research-development evidence**, not a substitute for results from a named physical machine with controlled power, affinity, and thermal state.

## Interpretation rules

1. Do not call an implementation universally “fastest” from these files.
2. Distinguish structural operation counts from measured latency.
3. Distinguish setup/allocation-inclusive call latency from planned execution-only FFT libraries.
4. Preserve raw samples and metadata when adding a result.
5. Treat tiny differences with uncertainty intervals spanning parity as statistically unresolved.
6. Update `SUMMARY.md` only with claims supported by checked-in data.

The full protocol is in [`../docs/EXPERIMENTS.md`](../docs/EXPERIMENTS.md).
