# Benchmark and accuracy results

This directory keeps **raw measurements and their provenance**, not just screenshots or headline numbers.

## Datasets

### `pr5-simd-kernel-baseline/`

The first formal machine-level codelet study. It holds the power-of-two radix-2 decomposition fixed while comparing:

- merged v3 `Radix2Plan`;
- scalar swap-list/stage-contiguous codelets;
- AVX2/FMA;
- AVX-512/FMA;
- experimental plan-time ISA selection;
- FFTW `ESTIMATE` and `MEASURE`.

The corpus contains **4,032 observations** across 3 randomized sessions and 6 sizes. It preserves execution and setup separately, records auto-selected ISA by session, and includes exact source/blob/binary/runtime provenance.

```bash
python3 tools/analyze_kernel.py results/pr5-simd-kernel-baseline/raw \
  --bootstrap 5000 --seed 20260812
```

### `pr4-fftw-baseline/`

The first controlled external production-library comparison. It compares fftlab reusable complex/real plans with FFTW 3.3.10 under matched setup/execution semantics and `ESTIMATE` / `MEASURE` planning. The corpus contains **3,456 observations** across 6 sizes.

```bash
python3 tools/analyze_vendor.py results/pr4-fftw-baseline/raw --seed 20260812
```

### `pr3-planning-real-baseline/`

The first formal dataset separating reusable setup from steady-state execution and studying a specialized real-input representation. It contains **2,790 timing observations** across six power-of-two sizes and five modes.

```bash
python3 tools/analyze_plan.py results/pr3-planning-real-baseline/raw --seed 20260812
```

### `pr2-research-baseline/`

The first formal expanded-algorithm dataset: raw randomized timing samples, structural complexity models, multi-family forward/backward numerical error, bootstrap uncertainty, and effect-size analysis.

### `baseline-linux-amd-epyc-gcc14.csv`

The v1 development baseline, retained for historical continuity.

## Interpretation rules

1. Do not call an implementation universally “fastest” from these files.
2. Distinguish structural operation counts from measured machine latency.
3. Distinguish one-time planning/setup from reusable execution cost.
4. Do not substitute complex-input FFT measurements for real-input workloads.
5. Production-library comparisons must normalize transform convention, inverse scaling, planner policy, allocation/workspace semantics, and thread count.
6. SIMD capability is not performance evidence: preserve scalar and every supported explicit-width result rather than reporting only an auto-selected path.
7. A wider ISA is not a winner unless the measured uncertainty supports that claim for the stated environment/size.
8. Auto-tuning cost and chosen policy are part of planning and must be reported; do not hide unsuccessful or unstable selections.
9. Preserve every raw sample and exact source/runtime provenance.
10. Treat uncertainty intervals spanning parity as statistically unresolved rather than forcing a winner.
11. Update `SUMMARY.md` only with claims supported by checked-in data.

The current formal baselines are from a containerized/virtualized environment. They support implementation-specific and methodology claims, not universal hardware rankings.

Protocols: [`../docs/EXPERIMENTS.md`](../docs/EXPERIMENTS.md), [`../docs/PLANNING_REAL.md`](../docs/PLANNING_REAL.md), [`../docs/VENDOR_BENCHMARKS.md`](../docs/VENDOR_BENCHMARKS.md), and [`../docs/SIMD_KERNELS.md`](../docs/SIMD_KERNELS.md).
