# Legacy research assets pending migration to EddyTodd/bench

The fft repository accumulated a useful empirical research layer while the algorithm implementations were still being designed. Version 1.0 separates that layer from the permanent installable library rather than deleting its evidence prematurely.

None of the files listed below are linked into `fftlab::fftlab` or installed by CMake.

## Benchmark / research C++ sources

These development-era sources should eventually move, be replaced by bench adapters, or be deleted after equivalent coverage exists in `EddyTodd/bench`:

- `src/research.cpp` — timing statistics, accuracy campaign glue, algorithm suites;
- `src/main.cpp` — historical research CLI;
- `src/plan_main.cpp` — plan timing CLI;
- `src/vendor_bench.cpp` — FFTW timing adapter;
- `src/kernel_bench.cpp` — SIMD timing/tuning study adapter;
- `src/arbitrary_bench.cpp` — arbitrary-plan/vendor timing adapter.

Several older implementation translation units (`src/common.cpp`, `src/power2.cpp`, `src/arbitrary.cpp`, `src/planned.cpp`, `src/arbitrary_plan.cpp`) remain development-history sources after the stable algorithms/plans were consolidated into the installed templated headers. They are not compiled by v1 and can be removed once historical benchmark adapters no longer depend on them.

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

Do not move or delete historical evidence until its exact source/raw-data relationship can be represented in `bench`. The v1 library completion milestone changes the build/API boundary now; evidence migration can then happen independently without destabilizing Fourier-transform code.
