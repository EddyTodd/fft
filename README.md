# fft

A dependency-free **C++23 Fourier-transform research laboratory** for studying FFT algorithms theoretically, numerically, and empirically.

The core algorithms are implemented from first principles rather than wrapped around FFTW, oneMKL, Accelerate/vDSP, cuFFT, or another tuned library. The public algorithm API lives in [`include/fftlab/fft.hpp`](include/fftlab/fft.hpp); reusable complex/real planning APIs live in [`include/fftlab/plan.hpp`](include/fftlab/plan.hpp); implementations are under [`src/`](src/). The surrounding repository defines the experimental methodology required to make performance and accuracy claims reproducible.

## Research questions

The project is built around questions that cannot be answered by Big-O notation alone:

- When does lower arithmetic complexity translate into lower wall-clock latency?
- What is the cost of recursion, allocation, bit reversal, autosorting, and extra workspace?
- How much latency is reusable setup/twiddle generation actually costing, and how quickly does planning amortize?
- How much faster should a real-input FFT be when Hermitian symmetry is exploited instead of ignored?
- For prime sizes, where does Rader beat Bluestein and where does it not?
- How strongly do factorization and transform size affect mixed-radix behavior?
- Which decomposition strategies accumulate the least floating-point error?
- How stable are rankings across sessions, compilers, and microarchitectures?
- When should an automatic dispatcher or planner prefer one implementation over another?

The general hypotheses and protocol are in [`docs/EXPERIMENTS.md`](docs/EXPERIMENTS.md). Mathematical context and structural operation models are in [`docs/THEORY.md`](docs/THEORY.md). Planning/execution semantics and the real-input reduction are specified separately in [`docs/PLANNING_REAL.md`](docs/PLANNING_REAL.md).

## Algorithm taxonomy

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

Forward and normalized inverse complex transforms are supported.

## Reusable plans and real-input transforms

`Radix2Plan` represents an N-point power-of-two transform whose reusable setup has already been performed. Construction stores:

- N bit-reversal indices;
- N/2 forward twiddle factors.

Its forward/inverse execution methods operate in place with **no dynamic allocation and no trigonometric setup**. This makes setup cost independently measurable instead of charging the same reusable work to every transform.

`RealRadix2Plan` specializes power-of-two real input. It packs the N real samples into one N/2-point planned complex FFT, applies O(N) recombination, and exposes exactly **N/2+1 nonredundant complex bins**. Its inverse consumes the packed half-spectrum and reconstructs N real samples.

The dedicated benchmark intentionally distinguishes:

- complex-plan construction;
- real-plan construction;
- legacy complex execution;
- planned complex execution;
- planned real execution.

See [`docs/PLANNING_REAL.md`](docs/PLANNING_REAL.md) for derivation, buffer semantics, benchmark contract, and amortization formulas.

## What makes this a research repository

The executables and tools provide more than algorithm implementations:

- deterministic input families: random, tones, impulse, alternating/cancellation-heavy, and high-dynamic-range;
- independent long-double O(N²) DFT reference;
- normalized forward and backward L1/L2/L∞ error measurements;
- forward→inverse round-trip error;
- structural complex-add/multiply and workspace models;
- explicit setup-vs-execution timing;
- reusable-plan and real-half-spectrum experiments;
- adaptive timing calibration and warmups;
- every raw timing sample retained on request;
- median, MAD, p05/p95, standard deviation, and deterministic bootstrap 95% median intervals;
- randomized multi-session experiment runners;
- pairwise bootstrap speedup intervals and effect sizes;
- checked-in raw results, environment metadata, hashes, source provenance, and analyses.

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

Complex algorithm suite:

```bash
./build/fft --self-test
```

The v2 complex-algorithm suite executes **14,396 checks** across prime, composite, arbitrary, and power-of-two sizes and all five signal families.

Planned/real suite:

```bash
./build/fft-plan --self-test
```

The dedicated v3 suite executes **1,233 checks** covering planned complex results, plan round trips, real half-spectrum equivalence, real round trips, and invalid domains. The exact v3 planned source has passed this suite under GCC 14.2, Clang 17, and GCC ASan/UBSan.

Sanitizers for the complete CMake project:

```bash
cmake -S . -B build-asan -DCMAKE_BUILD_TYPE=Debug -DFFT_SANITIZERS=ON
cmake --build build-asan -j
ctest --test-dir build-asan --output-on-failure
```

## Numerical accuracy

Inspect one complex algorithm/input:

```bash
./build/fft --verify --algorithm split-radix --size 1024 --signal random
./build/fft --verify --algorithm rader --size 1009 --signal dynamic-range
```

Generate a machine-readable study:

```bash
./build/fft --accuracy-suite --sizes 64,127,256,509 --csv > accuracy.csv
python3 tools/analyze_accuracy.py accuracy.csv
```

The methodology distinguishes **forward error** and **backward error** rather than relying on round-trip error alone. The reference is long double, not arbitrary precision; that limitation is explicit and remains a research frontier.

The real-plan correctness test independently verifies each returned half-spectrum bin against the corresponding nonredundant bin of a full complex FFT.

## Theoretical complexity

```bash
./build/fft --complexity --algorithm radix2-iterative --size 1024
./build/fft --complexity --algorithm split-radix --size 1024
./build/fft --complexity --algorithm rader --size 1009
```

These are structural complex-operation models, **not measured CPU FLOPs**. Twiddle generation, memory operations, indexing, vectorization, allocations, and microarchitectural effects are deliberately kept distinct from the mathematical model.

## Algorithm benchmark

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

## Planning and real-input benchmark

One size:

```bash
./build/fft-plan --benchmark --size 4096 --samples 31 --target-ms 5
```

A sweep:

```bash
./build/fft-plan --sweep --sizes 64,256,1024,4096,16384,65536 --csv
```

Preserve all five benchmark modes:

```bash
./build/fft-plan --sweep --sizes 64,256,1024 --samples 31 --raw-csv > plan-timings.csv
python3 tools/analyze_plan.py plan-timings.csv
```

Execution samples measure a forward+inverse pair and divide elapsed time by two. Buffers and plans are created before execution timing; plan construction is a separate mode. Mode order is randomized independently inside every sample.

Formal repeated experiment:

```bash
python3 tools/run_plan_experiment.py \
  --binary build/fft-plan \
  --out results/run-plan \
  --sizes 64,256,1024,4096,16384,65536 \
  --sessions 5 \
  --samples 51 \
  --target-ms 5 \
  --source-commit "$(git rev-parse HEAD)"

python3 tools/analyze_plan.py results/run-plan/timings.csv
```

The analyzer also accepts a directory containing plain `.csv` and gzip-compressed `.csv.gz` shards, which is how the checked-in v3 corpus is stored.

## General formal experiment

```bash
python3 tools/run_experiment.py \
  --binary build/fft \
  --sizes 64,127,256,509,1009,1024,4096 \
  --sessions 5 \
  --samples 51 \
  --target-ms 5
```

The runner records environment/protocol metadata, randomizes algorithm/size execution order separately per session, retains all raw timing observations, and captures bounded accuracy/complexity data.

See [`docs/EXPERIMENTS.md`](docs/EXPERIMENTS.md) before interpreting the output.

## Current empirical findings

### v3: planning and real input

The checked-in [`results/pr3-planning-real-baseline/`](results/pr3-planning-real-baseline/) dataset contains **2,790 raw observations**: 3 randomized sessions × 31 samples × 5 modes at 6 power-of-two sizes.

On the recorded container-exposed AMD EPYC 9V74 / GCC 14.2 environment:

- precomputed radix-2 execution is **1.295–1.593× faster** than the legacy equivalent path across N=64…65536; every bootstrap 95% interval excludes parity;
- complex-plan construction repays itself after roughly **1.75–3.95 transforms** according to median setup and execution savings;
- the specialized real-input path is **1.186× faster at N=64** and reaches roughly **1.74–1.75×** over planned complex execution at N=16384–65536;
- for N≥256, the measured extra real-plan setup cost is recovered within the first transform.

This establishes an important benchmark rule for subsequent milestones: **planning/setup, complex execution, and real execution are distinct workloads and must not be collapsed into one latency number**.

### v2: algorithm families and accuracy

The checked-in [`results/pr2-research-baseline/`](results/pr2-research-baseline/) dataset established, among other findings:

- iterative radix-2 is the practical power-of-two latency baseline in this implementation despite split-radix having a lower structural multiplication model;
- Rader is about **1.29× faster than Bluestein** at the recorded prime sizes N=509 and N=1009;
- at N=1024, direct DFT is roughly three orders of magnitude slower than iterative radix-2;
- split-radix and radix-4 have substantially lower measured forward L2 error than the iterative radix-2 implementation for the tested N=64 and N=256 input families;
- tiny wrapper/path differences are not called wins when uncertainty intervals include parity.

All current baselines are **descriptive findings for virtualized development environments**, not universal hardware rankings. See [`results/SUMMARY.md`](results/SUMMARY.md) for the evidence-backed headline history.

## Auto policy

The unplanned complex dispatcher remains intentionally inspectable:

- power of two → iterative radix-2;
- 2/3/5-smooth length → mixed-radix;
- prime N ≥ 17 → Rader;
- otherwise → Bluestein.

The reusable plan layer is currently narrower: it specializes power-of-two radix-2 complex and real transforms. Long term, planning should subsume fixed dispatch by allowing algorithm/codelet/back-end choices to be evaluated per hardware/compiler and workload while keeping plan construction separate from execution.

## Scope boundaries

The repository now covers single-threaded, double-precision, one-dimensional complex DFT algorithms plus reusable power-of-two complex and real radix-2 plans. It does **not yet** claim coverage of:

- arbitrary-length planned transforms;
- DCT/DST families;
- multidimensional FFTs;
- batched transforms;
- SIMD-specialized/generated codelets;
- multithreading;
- GPU kernels;
- persistent FFT “wisdom” or cross-process plan caches;
- arbitrary-precision reference arithmetic;
- controlled production-library comparisons with FFTW/pocketfft/oneMKL/vDSP/cuFFT/rocFFT;
- hardware-counter/cache/instruction analyses on named physical machines.

These are explicit future research milestones, not hidden omissions.

## Research integrity

A result should not enter the headline summary without raw data, environment metadata, exact source provenance, protocol/sample counts, uncertainty/effect-size analysis, and a scope statement. Benchmarks involving reusable libraries must additionally state whether planning, allocation, workspace preparation, and normalization are inside or outside the timed region.

See [`docs/EXPERIMENTS.md`](docs/EXPERIMENTS.md), [`docs/PLANNING_REAL.md`](docs/PLANNING_REAL.md), and [`CONTRIBUTING.md`](CONTRIBUTING.md).

## References

Primary literature and benchmark-methodology references are catalogued in [`docs/REFERENCES.md`](docs/REFERENCES.md).

## License

MIT. See [LICENSE](LICENSE).
