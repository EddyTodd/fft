# fft

A dependency-free **C++23 Fourier-transform research laboratory** for studying FFT algorithms theoretically, numerically, and empirically—and for comparing those first-principles implementations with optimized production libraries under controlled benchmark semantics.

The repository separates four questions that are often collapsed into one benchmark number:

1. **mathematical structure** — asymptotic and arithmetic complexity;
2. **numerical behavior** — forward/backward error across input families;
3. **implementation performance** — observed execution cost on a concrete machine;
4. **workload economics** — setup/planning cost plus repeated execution.

## Implemented algorithm families

| Algorithm | Domain | Research role |
|---|---|---|
| Direct DFT | any N | definition / quadratic baseline |
| Iterative radix-2 | powers of two | low-overhead baseline |
| Recursive radix-2 | powers of two | recursion/allocation study |
| Stockham radix-2 | powers of two | autosort / memory-traffic study |
| Radix-4 | powers of two | larger-butterfly study |
| Split-radix | powers of two | low arithmetic count vs actual latency |
| Mixed-radix | composite | factorization study |
| Rader | prime | prime DFT via cyclic convolution |
| Bluestein | any N | arbitrary-length convolution reduction |
| Auto | any N | inspectable empirical dispatch policy |

Reusable `Radix2Plan` and `RealRadix2Plan` APIs precompute permutation/twiddle state and expose allocation-free steady-state execution. The real plan returns only the `N/2+1` nonredundant bins implied by Hermitian symmetry.

## Research architecture

- [`include/fftlab/fft.hpp`](include/fftlab/fft.hpp) — algorithm API
- [`include/fftlab/plan.hpp`](include/fftlab/plan.hpp) — reusable complex/real plan API
- [`src/`](src/) — first-principles implementations and benchmark CLIs
- [`docs/THEORY.md`](docs/THEORY.md) — algorithm taxonomy and structural cost models
- [`docs/EXPERIMENTS.md`](docs/EXPERIMENTS.md) — experimental protocol and validity rules
- [`docs/PLANNING_REAL.md`](docs/PLANNING_REAL.md) — plan lifecycle and real-FFT derivation
- [`docs/VENDOR_BENCHMARKS.md`](docs/VENDOR_BENCHMARKS.md) — production-library normalization contract
- [`results/`](results/) — raw evidence, metadata, analyses, and headline history
- [`tools/`](tools/) — stdlib-only experiment runners and statistical analyzers

## Build and validation

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
```

For host-specific performance work:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DFFT_NATIVE=ON
cmake --build build -j
```

Sanitizers:

```bash
cmake -S . -B build-asan -DCMAKE_BUILD_TYPE=Debug -DFFT_SANITIZERS=ON
cmake --build build-asan -j
ctest --test-dir build-asan --output-on-failure
```

The complex-algorithm suite contains **14,396 checks**. The planned/real suite contains **1,233 dedicated checks**. The FFTW adapter additionally performs **495 direct frequency-bin cross-checks** against fftlab before benchmark claims are accepted.

## Numerical research

Five deterministic input families are available: random, multi-tone, impulse, alternating/cancellation-heavy, and high-dynamic-range.

```bash
./build/fft --verify --algorithm split-radix --size 1024 --signal random
./build/fft --accuracy-suite --sizes 64,127,256,509 --csv > accuracy.csv
python3 tools/analyze_accuracy.py accuracy.csv
```

Accuracy experiments use an independent long-double O(N²) DFT reference and report normalized forward/backward L1/L2/L∞ error plus round-trip error. Long double is explicitly treated as a limitation rather than arbitrary precision.

## Algorithm benchmarking

```bash
./build/fft --benchmark-suite --sizes 64,127,256,509,1009,1024
python3 tools/run_experiment.py --binary build/fft --sessions 5 --samples 51
```

The formal protocol preserves every raw timing observation, randomizes execution order, records exact source/environment provenance, and reports robust dispersion, bootstrap confidence intervals, and effect sizes.

## Planning and real-input benchmarking

```bash
./build/fft-plan --benchmark --size 4096 --samples 31 --target-ms 5
python3 tools/run_plan_experiment.py \
  --binary build/fft-plan \
  --out results/run-plan \
  --sizes 64,256,1024,4096,16384,65536 \
  --sessions 5 --samples 51 --target-ms 5
```

Planning and execution are separate response variables. Buffers and plans are created before execution timing; plan construction is measured independently. The checked-in v3 baseline showed reusable planning reducing steady-state radix-2 latency by **1.295–1.593×** and the specialized real path reaching roughly **1.74–1.75×** the planned complex throughput at larger tested sizes.

## FFTW production-library comparison

The `fft-vendor` executable dynamically loads a compatible FFTW double-precision runtime. FFTW is therefore an **optional runtime benchmark dependency**, not a dependency of the algorithm library itself.

```bash
./build/fft-vendor --info
./build/fft-vendor --self-test
./build/fft-vendor --raw-csv --size 4096 --samples 31 --setup-samples 1
```

Formal experiment:

```bash
python3 tools/run_vendor_experiment.py \
  --binary build/fft-vendor \
  --sizes 64,256,1024,4096,16384,65536 \
  --sessions 3 --samples 31 --setup-samples 1 \
  --source-commit "$(git rev-parse HEAD)" \
  --output results/run-vendor

python3 tools/analyze_vendor.py results/run-vendor/raw
```

The benchmark controls several easy-to-miss differences:

- complex is compared with complex; specialized real is compared with specialized real;
- both sides use persistent plans and preallocated execution buffers;
- setup/planning is timed separately;
- FFTW's unnormalized inverse receives its required `1/N` scaling **inside the timed execution path**;
- `FFTW_MEASURE` cold setup forgets wisdom before each measured plan pair;
- execution-mode and transform-size order are randomized;
- all raw samples are retained.

See [`docs/VENDOR_BENCHMARKS.md`](docs/VENDOR_BENCHMARKS.md).

## Current headline findings

### v4 — fftlab vs FFTW 3.3.10

The formal v4 corpus contains **3,456 raw observations** across 3 randomized sessions and six power-of-two sizes. On the recorded virtualized AMD EPYC 9V74 environment, compiled with GCC 14.2 and `-march=native`:

- FFTW `ESTIMATE` steady-state execution is **2.13–3.66× faster** than fftlab planned complex and **3.34–5.85× faster** than fftlab planned real across the recorded sizes;
- FFTW `MEASURE` reaches **3.77–4.68×** the fftlab complex execution speed and **4.16–6.59×** the fftlab real execution speed;
- that does **not** mean `MEASURE` is always the best workload choice: cold `MEASURE` plan creation costs tens of milliseconds to about two seconds in this run;
- the extra `MEASURE` planning effort needs roughly **19,668 to 2,126,344 transforms** to amortize versus `ESTIMATE`, depending on size and real/complex workload;
- at small N, fftlab's cheap plan can minimize total one-shot or short-run cost even though FFTW executes much faster; for complex N=64 the recorded `ESTIMATE` break-even is about **431 transforms**, falling to about **14 transforms at N=1024**;
- by complex N=16384 in this run, FFTW `ESTIMATE` is both cheaper to set up and faster to execute.

This is the repository's strongest demonstration so far that **“fastest FFT” is an incomplete question**. The answer depends on transform size, real vs complex data, planner policy, setup reuse count, library build, and hardware.

Full evidence: [`results/pr4-fftw-baseline/`](results/pr4-fftw-baseline/).

### Earlier milestones

- v3: planning and real-input specialization quantified separately from execution.
- v2: algorithm-family timing, structural complexity, numerical accuracy, Rader/Bluestein prime behavior, and statistical methodology.
- v1: dependency-free multi-algorithm baseline and benchmark harness.

See [`results/SUMMARY.md`](results/SUMMARY.md) for the evidence-backed history.

## Auto policy

The unplanned dispatcher remains intentionally inspectable:

- power of two → iterative radix-2;
- 2/3/5-smooth length → mixed-radix;
- prime N ≥ 17 → Rader;
- otherwise → Bluestein.

The long-term planner should supersede fixed dispatch where repeated execution justifies hardware/workload-specific selection.

## Scope and research frontier

The repository now covers one-dimensional, double-precision, single-threaded first-principles FFT algorithms, reusable power-of-two complex/real plans, and a controlled serial FFTW baseline. Important remaining work includes:

- Apple Accelerate/vDSP, Intel oneMKL, and pocketfft backends;
- matched alignment/vectorization studies;
- arbitrary-length reusable plans;
- `FFTW_PATIENT` / `EXHAUSTIVE` and wisdom persistence;
- multiple precisions and arbitrary-precision reference arithmetic;
- batched and multidimensional transforms;
- SIMD-specialized/generated codelets;
- multithreaded scaling;
- GPU backends;
- hardware performance counters and named physical-machine baselines across x86-64 and Arm.

## Research integrity

Headline results require raw data, environment metadata, exact source provenance, sample/session counts, uncertainty/effect-size analysis, and explicit timing semantics. Vendor comparisons additionally require normalization, planner policy, alignment/workspace, real-vs-complex, and threading semantics to be documented.

See [`docs/EXPERIMENTS.md`](docs/EXPERIMENTS.md), [`docs/VENDOR_BENCHMARKS.md`](docs/VENDOR_BENCHMARKS.md), and [`CONTRIBUTING.md`](CONTRIBUTING.md).

## License

MIT. See [LICENSE](LICENSE).
