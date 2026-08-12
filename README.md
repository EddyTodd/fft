# fft

A dependency-free **C++23 Fourier-transform research laboratory** for studying FFT algorithms theoretically and empirically.

This repository is intentionally not a wrapper around FFTW, oneMKL, Accelerate/vDSP, cuFFT, or another tuned library. The core algorithms are implemented from first principles in a small inspectable source tree, with the public API in [`include/fftlab/fft.hpp`](include/fftlab/fft.hpp) and implementations under [`src/`](src/), while the surrounding repository supplies the methodology required to make performance and numerical claims reproducible.

## Research questions

The project is built around questions that cannot be answered by Big-O notation alone:

- When does lower arithmetic complexity translate into lower wall-clock latency?
- What is the cost of recursion, allocation, bit reversal, autosorting, and extra workspace?
- For prime sizes, where does Rader beat Bluestein and where does it not?
- How strongly do factorization and transform size affect mixed-radix behavior?
- Which twiddle/decomposition strategies accumulate the least floating-point error?
- How stable are benchmark rankings across sessions, compilers, and microarchitectures?
- When should an automatic dispatcher prefer one algorithm over another?

The formal hypotheses and protocol are in [`docs/EXPERIMENTS.md`](docs/EXPERIMENTS.md). Mathematical context and operation models are in [`docs/THEORY.md`](docs/THEORY.md).

## Algorithms

| Algorithm | Supported N | Time | Key research role |
|---|---|---:|---|
| Direct DFT | any | O(N²) | definition baseline and independent cross-check |
| Iterative radix-2 Cooley–Tukey | powers of two | O(N log N) | low-overhead in-place baseline |
| Recursive radix-2 Cooley–Tukey | powers of two | O(N log N) | recursion/allocation experiment |
| Stockham radix-2 autosort | powers of two | O(N log N) | regular access vs extra workspace/traffic |
| Radix-4 Cooley–Tukey | powers of two | O(N log N) | fewer stages / larger butterflies |
| Split-radix | powers of two | O(N log N) | arithmetic-count vs practical-performance study |
| Mixed-radix Cooley–Tukey | composite; prime leaves allowed | typically O(N log N) | factorization experiments |
| Rader | prime | O(M log M) | prime DFT via cyclic convolution |
| Bluestein / chirp-z | any | O(M log M) | robust arbitrary-length convolution reduction |
| Auto dispatcher | any | varies | inspectable, empirically revisable policy |

Forward and normalized inverse transforms are supported.

## What makes this a research repository

The executable contains more than algorithm implementations:

- deterministic input families: random, tones, impulse, alternating/cancellation-heavy, and high-dynamic-range;
- independent long-double O(N²) DFT reference;
- normalized forward and backward L1/L2/L∞ error measurements;
- forward→inverse round-trip error;
- structural complex-add/multiply and workspace models;
- adaptive timing calibration and warmups;
- every raw timing sample on request;
- median, MAD, p05/p95, standard deviation, and deterministic bootstrap 95% median CI;
- a randomized multi-session experiment runner;
- pairwise bootstrap speedup intervals and common-language effect sizes;
- checked-in raw results, environment metadata, and analyses.

## Build

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
```

Requires CMake 3.20+ and a C++23 compiler.

For a host-specific research build:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DFFT_NATIVE=ON -DFFT_LTO=ON
cmake --build build -j
```

`FFT_NATIVE=ON` improves relevance to one machine but reduces portability of the resulting binary and benchmark.

## Correctness

```bash
./build/fft --self-test
```

The current integrated suite executes **14,396 checks** across prime, composite, arbitrary, and power-of-two sizes and all five signal families. It cross-checks supported algorithms against the direct DFT and verifies inverse behavior and domain rejection.

Sanitizers:

```bash
cmake -S . -B build-asan -DCMAKE_BUILD_TYPE=Debug -DFFT_SANITIZERS=ON
cmake --build build-asan -j
ctest --test-dir build-asan --output-on-failure
```

## Numerical accuracy

Inspect one algorithm/input:

```bash
./build/fft --verify --algorithm split-radix --size 1024 --signal random
./build/fft --verify --algorithm rader --size 1009 --signal dynamic-range
```

Generate a machine-readable study:

```bash
./build/fft --accuracy-suite --sizes 64,127,256,509 --csv > accuracy.csv
python3 tools/analyze_accuracy.py accuracy.csv
```

The methodology follows the important distinction between **forward error** and **backward error** rather than relying on round-trip error alone. The reference is long double, not arbitrary precision; that limitation is explicit and is part of the research roadmap.

## Theoretical complexity

```bash
./build/fft --complexity --algorithm radix2-iterative --size 1024
./build/fft --complexity --algorithm split-radix --size 1024
./build/fft --complexity --algorithm rader --size 1009
```

These are structural complex-operation models, **not measured CPU FLOPs**. Twiddle generation, memory operations, indexing, vectorization, allocations, and microarchitectural effects are deliberately kept distinct from the mathematical model.

## Quick benchmark

```bash
./build/fft --benchmark --algorithm radix2-iterative --size 4096
./build/fft --benchmark --algorithm rader --size 1009 --samples 51 --target-ms 10
./build/fft --benchmark-suite --sizes 64,127,256,509,1009,1024
```

Preserve raw samples:

```bash
./build/fft --benchmark-suite --sizes 64,127,256 --samples 51 --raw-csv > timings.csv
python3 tools/analyze.py timings.csv
```

The quick suite runs algorithms sequentially. Use the formal experiment runner for publishable comparisons.

## Formal experiment

```bash
python3 tools/run_experiment.py \
  --binary build/fft \
  --sizes 64,127,256,509,1009,1024,4096 \
  --sessions 5 \
  --samples 51 \
  --target-ms 5
```

The runner:

1. records environment/protocol metadata;
2. randomizes `(N, algorithm)` execution order separately per session;
3. retains all raw timing observations;
4. captures accuracy data for bounded sizes;
5. captures theoretical complexity rows.

Then:

```bash
python3 tools/analyze.py results/run-*/timings.csv > timing-analysis.md
python3 tools/analyze_accuracy.py results/run-*/accuracy.csv > accuracy-analysis.md
```

See [`docs/EXPERIMENTS.md`](docs/EXPERIMENTS.md) before interpreting the output.

## Current empirical findings

The checked-in v2 research baseline contains 3 randomized sessions × 31 samples per algorithm/size on the container-exposed AMD EPYC 9V74 environment. Important observations include:

- iterative radix-2 remains the practical power-of-two latency baseline in this implementation despite split-radix having a more attractive structural arithmetic model;
- at prime N=509, Rader’s median latency is about 1.31× lower than Bluestein’s in the recorded run;
- at N=1024, direct DFT is roughly three orders of magnitude slower than iterative radix-2;
- split-radix and radix-4 have substantially lower measured forward L2 error than the iterative radix-2 implementation at the tested N=64 and N=256 inputs;
- the `auto` wrapper and its selected implementation are often statistically indistinguishable at small speedup margins, as expected from identical transform paths plus measurement noise.

These are **descriptive findings for one virtualized environment**, not universal rankings. Raw data and uncertainty analysis live under [`results/pr2-research-baseline/`](results/pr2-research-baseline/).

## Auto policy

The current policy is intentionally inspectable:

- power of two → iterative radix-2;
- 2/3/5-smooth length → mixed-radix;
- prime N ≥ 17 → Rader;
- otherwise → Bluestein.

This is a research baseline, not a claim of global optimality. The correct long-term direction is an FFTW-style planning layer that benchmarks/learns crossover decisions per hardware/compiler and separates planning cost from repeated execution cost.

## Scope boundaries

The repository currently studies **single-threaded, double-precision, one-dimensional complex DFTs**. It does not yet claim coverage of:

- real-input FFT specializations;
- DCT/DST families;
- multidimensional FFTs;
- batched transforms;
- SIMD-specialized kernels;
- multithreading;
- GPU kernels;
- plan reuse/precomputed twiddle tables;
- arbitrary-precision reference arithmetic;
- production-library comparisons with FFTW/oneMKL/vDSP/cuFFT/rocFFT.

Those are explicit future research milestones, not hidden omissions.

## Research integrity

A result should not enter the headline summary without raw data, environment metadata, exact source provenance, protocol/sample counts, uncertainty/effect-size analysis, and a scope statement. See [`docs/EXPERIMENTS.md`](docs/EXPERIMENTS.md) and [`CONTRIBUTING.md`](CONTRIBUTING.md).

## References

Primary literature and benchmark-methodology references are catalogued in [`docs/REFERENCES.md`](docs/REFERENCES.md).

## License

MIT. See [LICENSE](LICENSE).
