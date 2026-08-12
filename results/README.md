# Benchmark and accuracy results

This directory keeps **raw measurements and their provenance**, not just screenshots or headline numbers.

## Datasets

### `pr6-arbitrary-plan-baseline/`

The first formal reusable-plan study for non-power-of-two complex FFTs. It compares legacy and planned Bluestein/Rader reductions with persistent FFTW `ESTIMATE` / `MEASURE` plans at eight prime sizes.

The evidence represents **4,560 raw observations** across 3 randomized sessions. It separates setup from execution, records convolution length and Rader's direct-cyclic special case, and retains exact source/blob/binary/runtime provenance.

Because this PR was written through a text-only connector, each canonical gzip CSV is checked in **losslessly as base64 parts** under `raw/`. `tools/analyze_arbitrary.py` detects and reconstructs those parts in memory. `metadata.json` records the SHA-256 of each original gzip stream, so the transport representation can be verified independently.

```bash
python3 tools/analyze_arbitrary.py results/pr6-arbitrary-plan-baseline/raw \
  --bootstrap 5000 --seed 20260812
```

### `pr5-simd-kernel-baseline/`

The first formal machine-level codelet study. It holds the power-of-two radix-2 decomposition fixed while comparing the merged v3 plan, scalar layout, AVX2/FMA, AVX-512/FMA, experimental plan-time ISA selection, and FFTW. The corpus contains **4,032 observations**.

```bash
python3 tools/analyze_kernel.py results/pr5-simd-kernel-baseline/raw \
  --bootstrap 5000 --seed 20260812
```

### `pr4-fftw-baseline/`

The first controlled external production-library comparison. It compares fftlab reusable complex/real plans with FFTW 3.3.10 under matched setup/execution semantics and `ESTIMATE` / `MEASURE` planning. The corpus contains **3,456 observations**.

### `pr3-planning-real-baseline/`

The first formal dataset separating reusable setup from steady-state execution and studying a specialized real-input representation. It contains **2,790 observations**.

### `pr2-research-baseline/`

The first formal expanded-algorithm dataset: randomized timing samples, structural complexity models, multi-family forward/backward numerical error, bootstrap uncertainty, and effect-size analysis.

### `baseline-linux-amd-epyc-gcc14.csv`

The v1 development baseline, retained for historical continuity.

## Interpretation rules

1. Do not call an implementation universally “fastest” from these files.
2. Distinguish structural operation/convolution sizes from measured machine latency.
3. Distinguish one-time planning/setup from reusable execution cost.
4. Compare legacy setup-inclusive APIs only as APIs; use planned-vs-planned data for reusable algorithm ranking.
5. Do not substitute complex-input measurements for specialized real-input workloads.
6. Production-library comparisons must normalize transform convention, inverse scaling, planner policy, allocation/workspace semantics, and thread count.
7. SIMD capability is not performance evidence; preserve scalar and every supported explicit-width result.
8. A wider ISA is not a winner unless measured uncertainty supports the claim for the stated environment/size.
9. Auto-tuning cost and selected policy are part of planning; preserve unsuccessful or unstable selections.
10. Prime dispatch must remain auditable: preserve explicit Bluestein and Rader results rather than reporting only the automatic choice.
11. Preserve every raw sample and exact source/runtime provenance, including any reversible transport encoding used to store it.
12. Treat intervals spanning parity as statistically unresolved rather than forcing a winner.
13. Update `SUMMARY.md` only with claims supported by checked-in evidence.

The current formal baselines are from a containerized/virtualized environment. They support implementation-specific and methodology claims, not universal hardware rankings.

Protocols: [`../docs/EXPERIMENTS.md`](../docs/EXPERIMENTS.md), [`../docs/PLANNING_REAL.md`](../docs/PLANNING_REAL.md), [`../docs/ARBITRARY_PLANS.md`](../docs/ARBITRARY_PLANS.md), [`../docs/VENDOR_BENCHMARKS.md`](../docs/VENDOR_BENCHMARKS.md), and [`../docs/SIMD_KERNELS.md`](../docs/SIMD_KERNELS.md).
