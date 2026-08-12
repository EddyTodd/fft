# fft

A dependency-free **C++23 Fourier-transform research laboratory** for studying FFT algorithms theoretically, numerically, and empirically—and comparing first-principles implementations with optimized production libraries under controlled benchmark semantics.

The project deliberately separates questions that are often collapsed into one timing number:

1. **mathematical structure** — asymptotic complexity and arithmetic/decomposition choices;
2. **numerical behavior** — forward/backward error across representative and adversarial inputs;
3. **implementation performance** — latency, data movement, layout, code generation, and SIMD width;
4. **planning economics** — one-time setup versus reusable steady-state execution;
5. **workload specialization** — real/complex data, prime/composite structure, planner policy, and hardware-specific choices.

## Algorithm and implementation taxonomy

| Family / implementation | Domain | Research role |
|---|---|---|
| Direct DFT | any N | definition and quadratic baseline |
| Iterative / recursive radix-2 | powers of two | iterative overhead and recursion/allocation studies |
| Stockham radix-2 | powers of two | autosort and memory-traffic study |
| Radix-4 / split-radix | powers of two | butterfly/arithmetic-count studies |
| Mixed-radix | composite | factorization study |
| Rader | prime | prime DFT via cyclic convolution |
| Bluestein | any N | arbitrary-length convolution reduction |
| `Radix2Plan` | powers of two | reusable permutation/twiddle state |
| `RealRadix2Plan` | powers of two | packed `N/2+1` real spectrum |
| `BluesteinPlan` | any N >= 1 | reusable chirp + convolution-kernel spectrum |
| `RaderPlan` | prime N >= 3 | reusable prime permutations + convolution kernel |
| `ArbitraryPlan` | any N >= 1 | structural reusable planner across radix-2/Rader/Bluestein |
| `KernelRadix2Plan` scalar / AVX2 / AVX-512 | x86 power-of-two | matched scalar/SIMD codelet research |
| `KernelRadix2Plan::Auto` | supported power-of-two | experimental plan-time empirical ISA selection |
| FFTW `ESTIMATE` / `MEASURE` | runtime dependent | controlled production-library baseline |

The core algorithms are implemented from first principles. FFTW is dynamically loaded only by benchmark executables and is **not** a dependency of the fftlab library.

## Repository map

- [`include/fftlab/fft.hpp`](include/fftlab/fft.hpp) — algorithm API
- [`include/fftlab/plan.hpp`](include/fftlab/plan.hpp) — reusable power-of-two complex/real plans
- [`include/fftlab/arbitrary_plan.hpp`](include/fftlab/arbitrary_plan.hpp) — reusable Bluestein/Rader/general plans
- [`include/fftlab/kernel.hpp`](include/fftlab/kernel.hpp) — runtime-selectable radix-2 codelet plan
- [`docs/EXPERIMENTS.md`](docs/EXPERIMENTS.md) — general empirical protocol
- [`docs/PLANNING_REAL.md`](docs/PLANNING_REAL.md) — planning lifecycle and real FFT derivation
- [`docs/ARBITRARY_PLANS.md`](docs/ARBITRARY_PLANS.md) — arbitrary-length plan derivation, memory model, and prime protocol
- [`docs/VENDOR_BENCHMARKS.md`](docs/VENDOR_BENCHMARKS.md) — external-library normalization contract
- [`docs/SIMD_KERNELS.md`](docs/SIMD_KERNELS.md) — codelet/ISA design and validity rules
- [`results/`](results/) — raw evidence, metadata, analyses, and headline history
- [`tools/`](tools/) — standard-library experiment runners and analyzers

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

Validation depth now includes:

- historical complex-algorithm suite: **14,396 checks**;
- planned/real suite: **1,233 checks**;
- v5 exact kernel source: **11,021 SIMD/codelet checks + 328 FFTW cross-checks**;
- v6 final arbitrary-plan policy: **5,574 arbitrary-plan checks + 653 FFTW prime-bin cross-checks**.

The v5 and v6 specialized suites passed optimized GCC 14.2, optimized Clang 17, and GCC ASan/UBSan in the recorded development environment.

## Core research workflows

Algorithm families and accuracy:

```bash
./build/fft --verify --algorithm split-radix --size 1024 --signal random
./build/fft --benchmark-suite --sizes 64,127,256,509,1009,1024
python3 tools/run_experiment.py --binary build/fft --sessions 5 --samples 51
```

Reusable power-of-two planning and real input:

```bash
./build/fft-plan --benchmark --size 4096 --samples 31 --target-ms 5
python3 tools/run_plan_experiment.py \
  --binary build/fft-plan --out results/run-plan \
  --sizes 64,256,1024,4096,16384,65536 \
  --sessions 5 --samples 51 --target-ms 5
```

Production-library comparison:

```bash
./build/fft-vendor --info
./build/fft-vendor --self-test
```

SIMD/codelet study:

```bash
./build/fft-kernel --info
./build/fft-kernel --self-test
python3 tools/run_kernel_experiment.py \
  --binary build/fft-kernel --out results/run-kernel \
  --sizes 64,256,1024,4096,16384,65536 \
  --sessions 3 --samples 31 --setup-samples 1 --target-ms 2 \
  --source-commit "$(git rev-parse HEAD)"
```

### Arbitrary-length reusable plans

`fft-arbitrary` compares setup-inclusive historical APIs with reusable Bluestein/Rader plans and FFTW under one prime-length harness:

```bash
./build/fft-arbitrary --info
./build/fft-arbitrary --self-test
./build/fft-arbitrary --raw-csv --size 509 --samples 31 --setup-samples 1
```

Formal study:

```bash
python3 tools/run_arbitrary_experiment.py \
  --binary build/fft-arbitrary \
  --out results/run-arbitrary \
  --sizes 17,31,61,127,257,509,1009,4093 \
  --sessions 3 --samples 31 --setup-samples 1 --target-ms 2 \
  --source-commit "$(git rev-parse HEAD)"

python3 tools/analyze_arbitrary.py results/run-arbitrary/raw \
  --bootstrap 5000 --seed 20260812
```

`BluesteinPlan` precomputes chirps, an FFT-domain convolution kernel, and an M-point radix-2 plan. `RaderPlan` precomputes prime permutations and its convolution-kernel spectrum. When `p-1` is a power of two, Rader performs the cyclic convolution **directly at length `p-1`** rather than zero-padding it to a linear convolution.

`ArbitraryPlan::Auto` uses a structural policy rather than the old blanket `prime -> Rader` rule:

- power of two -> radix-2;
- prime where planned Rader has a **strictly shorter convolution** than Bluestein -> Rader;
- otherwise -> Bluestein.

Explicit Rader and Bluestein policies remain available so the decision can be audited and retested.

## Current headline findings

### v6 — arbitrary-length planning and prime crossover

The checked-in [`results/pr6-arbitrary-plan-baseline/`](results/pr6-arbitrary-plan-baseline/) evidence represents **4,560 raw observations** across 3 randomized sessions and eight prime sizes.

On the recorded virtualized AMD EPYC / GCC 14.2 environment:

- reusable planning improves Bluestein by roughly **3.09–4.47×** versus the historical setup-inclusive API;
- reusable planning improves Rader by roughly **2.35–7.61×**;
- plan construction generally amortizes in about **0.8–3.3 transforms**;
- at **N=17**, Rader's direct 16-point cyclic convolution is **4.34× faster** than planned Bluestein's 64-point convolution;
- at **N=257**, Rader's direct 256-point cyclic convolution is **4.33× faster** than planned Bluestein's 1024-point convolution;
- when both reductions use the same power-of-two convolution length, the differences become small: Rader wins N=31 by ~1.8%, Bluestein wins N=61, 127, 509, and 1009 by roughly 0.8–1.6%, and N=4093 is statistically unresolved;
- the best first-principles planned reduction is about **1.28–3.65× slower** than FFTW `MEASURE` across the formal prime matrix.

This materially revises the earlier v2 observation that Rader beat Bluestein at N=509 and 1009: those measurements used the legacy setup-inclusive APIs. Once reusable planning is normalized, **Bluestein is slightly faster at both sizes in the v6 corpus**. That is a central research result: lifecycle semantics can change an apparent algorithm ranking.

Full derivation and buffer/memory contracts: [`docs/ARBITRARY_PLANS.md`](docs/ARBITRARY_PLANS.md).

### v5 — SIMD radix-2 codelets

The v5 corpus contains **4,032 observations**. The best explicit SIMD codelet is **1.747–2.104× faster** than the merged v3 plan and closes roughly **57.9–68.2%** of its latency gap to FFTW `MEASURE`. AVX2 wins the recorded small/medium sizes while AVX-512 wins the larger sizes; the experimental auto-tuner is intentionally documented as imperfect.

### v4 — FFTW production baseline

The v4 corpus contains **3,456 observations**. FFTW `MEASURE` is much faster in steady-state execution, but cold planning can cost tens of milliseconds to seconds and require tens of thousands to millions of transforms to amortize over `ESTIMATE`.

### Earlier milestones

- v3 — reusable planning and specialized real-input transforms;
- v2 — broad algorithm timing, structural complexity, numerical accuracy, and initial Rader/Bluestein observations;
- v1 — dependency-free multi-algorithm baseline.

See [`results/SUMMARY.md`](results/SUMMARY.md) for the evidence-backed history.

## Research frontier

Important next work includes:

- reusable **mixed-radix convolution plans** so Rader is not forced through power-of-two zero padding when `p-1` factors well;
- generated small-N/mixed-radix codelets and schedule selection;
- vectorized real-input recombination and arbitrary-length SIMD kernels;
- cache-blocked, four-step, six-step, and cache-oblivious large transforms;
- Apple Accelerate/vDSP, Intel oneMKL, and pocketfft backends;
- Arm NEON/SVE plus named physical x86-64/Arm baselines;
- hardware performance counters;
- persisted/statistically robust tuning wisdom;
- `FFTW_PATIENT` / `EXHAUSTIVE` and wisdom persistence;
- multiple precisions and arbitrary-precision reference arithmetic;
- batched, multidimensional, multithreaded, and GPU transforms.

## Research integrity

Headline results require raw data, exact source/binary/runtime provenance, sample/session counts, uncertainty analysis, and explicit timing semantics. Planned comparisons must state what is persistent and what scratch is caller-owned. Vendor comparisons additionally require planner, normalization, allocation/alignment, representation, and threading semantics. Architecture claims require explicit ISA controls and runtime capability gates.

See [`docs/EXPERIMENTS.md`](docs/EXPERIMENTS.md), [`docs/ARBITRARY_PLANS.md`](docs/ARBITRARY_PLANS.md), [`docs/SIMD_KERNELS.md`](docs/SIMD_KERNELS.md), [`docs/VENDOR_BENCHMARKS.md`](docs/VENDOR_BENCHMARKS.md), and [`CONTRIBUTING.md`](CONTRIBUTING.md).

## License

MIT. See [LICENSE](LICENSE).
