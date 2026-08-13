# Legacy research assets pending migration to EddyTodd/bench

The FFT repository accumulated a useful empirical research layer while the algorithm implementations were still being designed. Version 1.0 separates that layer physically as well as logically from the permanent installable library rather than deleting its evidence prematurely.

None of the files listed below are linked into `fftlab::fftlab` or installed by CMake.

## Benchmark / research C++ sources

Historical frontends now live under `research/apps/`:

- `main.cpp` — historical research CLI;
- `plan_main.cpp` — plan timing CLI;
- `vendor_bench.cpp` — FFTW timing adapter;
- `kernel_bench.cpp` — SIMD timing/tuning adapter;
- `arbitrary_bench.cpp` — arbitrary-plan/vendor timing adapter.

Development-era implementation/support sources now live under `research/legacy-src/`:

- `research.cpp` — timing statistics, accuracy campaign glue, algorithm suites;
- `common.cpp`;
- `power2.cpp`;
- `arbitrary.cpp`;
- `planned.cpp`;
- `arbitrary_plan.cpp`;
- `internal.hpp`.

These remain byte-identical historical sources. They are retained only until corresponding benchmark/evidence migration is complete.

The permanent `src/` directory now contains only `kernel.cpp`, the non-header implementation compiled by v1. The stable algorithms and plans otherwise live in installed templated headers.

## Python campaign and analysis tools

Move to `EddyTodd/bench` or replace with generic bench equivalents:

- `tools/run_experiment.py`
- `tools/run_plan_experiment.py`
- `tools/run_vendor_experiment.py`
- `tools/run_kernel_experiment.py`
- `tools/run_arbitrary_experiment.py`
- `tools/analyze.py`
- `tools/analyze_accuracy.py`
- `tools/analyze_plan.py`
- `tools/analyze_vendor.py`
- `tools/analyze_kernel.py`
- `tools/analyze_arbitrary.py`

## Results/evidence

Historical result corpora under `results/` should migrate as immutable benchmark evidence or be archived once `bench` owns their schemas:

- `results/baseline-linux-amd-epyc-gcc14.csv`
- `results/pr2-research-baseline/`
- `results/pr3-planning-real-baseline/`
- `results/pr4-fftw-baseline/`
- `results/pr5-simd-kernel-baseline/`
- `results/pr6-arbitrary-plan-baseline/`
- `results/SUMMARY.md`

## Research-methodology documents

These are valuable context but are not permanent core API specifications:

- `docs/EXPERIMENTS.md`
- `docs/VENDOR_BENCHMARKS.md`
- result-specific portions of `docs/SIMD_KERNELS.md`, `docs/PLANNING_REAL.md`, and `docs/ARBITRARY_PLANS.md`.

During migration, preserve any FFT-specific mathematical explanation in this repository and move generic sampling/statistics/provenance/campaign material to `bench`.

## Migration rule

Do not delete historical evidence until its exact source/raw-data relationship can be represented in `bench`. Physical relocation inside this repository preserves the same blobs and Git history; it is source-tree sanitation, not evidence deletion.
