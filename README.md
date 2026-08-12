# fft

A dependency-free **C++23 Fourier-transform research laboratory** for studying FFT algorithms theoretically, numerically, and empirically—and comparing first-principles implementations with optimized production libraries under controlled benchmark semantics.

The project separates questions that are often collapsed into a single timing number:

1. **mathematical structure** — asymptotic complexity and structural arithmetic;
2. **numerical behavior** — forward/backward error across adversarial and representative inputs;
3. **implementation performance** — latency, data movement, code generation, SIMD width, and memory layout;
4. **planning economics** — one-time setup versus reusable steady-state execution;
5. **workload specialization** — complex versus real input, planner policy, and hardware-specific choices.

## Algorithm and implementation taxonomy

| Family / implementation | Domain | Research role |
|---|---|---|
| Direct DFT | any N | definition and quadratic baseline |
| Iterative radix-2 | powers of two | low-overhead algorithm baseline |
| Recursive radix-2 | powers of two | recursion/allocation study |
| Stockham radix-2 | powers of two | autosort and memory-traffic study |
| Radix-4 | powers of two | larger-butterfly study |
| Split-radix | powers of two | arithmetic count versus actual latency |
| Mixed-radix | composite | factorization study |
| Rader | prime | prime DFT via cyclic convolution |
| Bluestein | any N | arbitrary-length convolution reduction |
| `Radix2Plan` | powers of two | reusable permutation/twiddle state |
| `RealRadix2Plan` | powers of two | packed `N/2+1` real spectrum |
| `KernelRadix2Plan` scalar | powers of two | swap-list/stage-local codelet baseline |
| `KernelRadix2Plan` AVX2/FMA | x86 power-of-two | two-complex SIMD codelet |
| `KernelRadix2Plan` AVX-512/FMA | x86 power-of-two | four-complex SIMD codelet |
| `KernelRadix2Plan` Auto | supported power-of-two | experimental plan-time empirical ISA selection |
| FFTW `ESTIMATE` / `MEASURE` | runtime dependent | controlled production-library baseline |

The core algorithms remain implemented from first principles. FFTW is loaded only by benchmark executables at runtime and is not a dependency of the fftlab algorithm library.

## Repository map

- [`include/fftlab/fft.hpp`](include/fftlab/fft.hpp) — algorithm API
- [`include/fftlab/plan.hpp`](include/fftlab/plan.hpp) — reusable complex/real plans
- [`include/fftlab/kernel.hpp`](include/fftlab/kernel.hpp) — runtime-selectable radix-2 codelet plan
- [`src/`](src/) — implementations and benchmark CLIs
- [`docs/THEORY.md`](docs/THEORY.md) — algorithm taxonomy and structural models
- [`docs/EXPERIMENTS.md`](docs/EXPERIMENTS.md) — general empirical protocol
- [`docs/PLANNING_REAL.md`](docs/PLANNING_REAL.md) — planning lifecycle and real FFT derivation
- [`docs/VENDOR_BENCHMARKS.md`](docs/VENDOR_BENCHMARKS.md) — external-library comparison contract
- [`docs/SIMD_KERNELS.md`](docs/SIMD_KERNELS.md) — SIMD/codelet design, tuning, memory model, and validity rules
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

The historical complex-algorithm suite contains **14,396 checks** and the planned/real suite **1,233 checks**. The v5 kernel milestone adds **11,021 SIMD/codelet checks** plus **328 independent FFTW frequency-bin cross-checks** for the exact benchmark source. The kernel checks passed under GCC 14.2, Clang 17, and GCC ASan/UBSan.

## Numerical research

```bash
./build/fft --verify --algorithm split-radix --size 1024 --signal random
./build/fft --accuracy-suite --sizes 64,127,256,509 --csv > accuracy.csv
python3 tools/analyze_accuracy.py accuracy.csv
```

Five deterministic signal families are used. Accuracy experiments use an independent long-double `O(N²)` DFT reference and report normalized forward/backward L1/L2/L∞ error plus round-trip error. Long double is explicitly treated as a limitation rather than arbitrary precision.

## Algorithm-family benchmarking

```bash
./build/fft --benchmark-suite --sizes 64,127,256,509,1009,1024
python3 tools/run_experiment.py --binary build/fft --sessions 5 --samples 51
```

The protocol preserves every raw timing observation, randomizes execution order, records source/environment provenance, and reports robust dispersion, bootstrap intervals, and effect sizes.

## Planning and real-input benchmarking

```bash
./build/fft-plan --benchmark --size 4096 --samples 31 --target-ms 5
python3 tools/run_plan_experiment.py \
  --binary build/fft-plan \
  --out results/run-plan \
  --sizes 64,256,1024,4096,16384,65536 \
  --sessions 5 --samples 51 --target-ms 5
```

Planning and execution are separate response variables. The v3 baseline showed reusable planning reducing steady-state radix-2 latency by **1.295–1.593×** and the specialized real path reaching roughly **1.74–1.75×** the planned complex throughput at the larger recorded sizes.

## Production-library benchmarking

`fft-vendor` dynamically loads a compatible FFTW double-precision runtime:

```bash
./build/fft-vendor --info
./build/fft-vendor --self-test
./build/fft-vendor --raw-csv --size 4096 --samples 31 --setup-samples 1
```

The v4 benchmark compares persistent plans with matched normalization, real/complex representation, allocation semantics, one thread, randomized order, separately measured planning, and retained raw samples. See [`docs/VENDOR_BENCHMARKS.md`](docs/VENDOR_BENCHMARKS.md).

## SIMD/codelet benchmark

`fft-kernel` holds the radix-2 decomposition fixed while comparing the merged reusable plan, a scalar codelet layout, AVX2/FMA, AVX-512/FMA, experimental auto selection, and FFTW in one harness.

```bash
./build/fft-kernel --info
./build/fft-kernel --self-test
./build/fft-kernel --raw-csv --size 4096 --samples 31 --setup-samples 1
```

Formal multi-session experiment:

```bash
python3 tools/run_kernel_experiment.py \
  --binary build/fft-kernel \
  --out results/run-kernel \
  --sizes 64,256,1024,4096,16384,65536 \
  --sessions 3 --samples 31 --setup-samples 1 \
  --target-ms 2 \
  --source-commit "$(git rev-parse HEAD)"

python3 tools/analyze_kernel.py results/run-kernel/raw \
  --bootstrap 5000 --seed 20260812
```

The explicit SIMD functions use function-level ISA targets and runtime capability checks, so the portable library build does not globally require AVX2 or AVX-512. Unsupported explicit kernels are rejected rather than executed speculatively.

## Current headline findings

### v5 — scalar, AVX2, AVX-512, and empirical ISA selection

The checked-in [`results/pr5-simd-kernel-baseline/`](results/pr5-simd-kernel-baseline/) corpus contains **4,032 raw observations** from 3 randomized sessions at N=64…65536.

On the recorded virtualized AMD EPYC 9V74 / GCC 14.2 environment:

- the best explicit SIMD codelet is **1.747–2.104× faster** than the merged v3 reusable plan;
- that one optimization layer closes roughly **57.9–68.2%** of the latency gap from the v3 plan to FFTW `MEASURE`;
- FFTW `MEASURE` still remains approximately **1.97–2.35× faster** than the best local SIMD codelet;
- AVX2 is significantly faster than AVX-512 at N=64, 256, and 1024;
- AVX-512 is significantly faster at N=4096, 16384, and 65536;
- every AVX2-versus-AVX-512 bootstrap interval in the formal matrix excludes parity;
- explicit SIMD setup amortizes quickly—about **20 transforms at N=64**, falling below one transform at N≥16384 in the recorded medians;
- the experimental `Auto` tuner costs roughly **26–99 ms** to construct and does **not** always choose the pooled explicit-kernel winner at intermediate sizes.

The last point is intentionally preserved as a result rather than tuned away. A supported or wider ISA is not proof that it is fastest, and a lightweight noisy tuner is not automatically a reliable production dispatch policy.

Full derivation, memory tradeoffs, capability rules, and threats to validity: [`docs/SIMD_KERNELS.md`](docs/SIMD_KERNELS.md).

### v4 — fftlab versus FFTW 3.3.10

The v4 corpus contains **3,456 raw observations**. FFTW `MEASURE` reached roughly **3.77–4.68×** the fftlab planned-complex execution speed and **4.16–6.59×** the planned-real speed, but cold `MEASURE` planning required tens of milliseconds to about two seconds. Depending on N/workload, roughly **19,668 to 2,126,344 transforms** were needed to amortize `MEASURE` over `ESTIMATE`.

That result established that steady-state throughput and total workload cost can have different winners.

### Earlier milestones

- v3 — reusable planning and real-input specialization quantified separately from execution;
- v2 — algorithm-family timing, structural complexity, numerical accuracy, Rader/Bluestein behavior, and statistical methodology;
- v1 — dependency-free multi-algorithm baseline and benchmark harness.

See [`results/SUMMARY.md`](results/SUMMARY.md) for the evidence-backed history.

## Dispatch policy

Two different dispatch questions are deliberately kept separate:

- the **unplanned algorithm dispatcher** chooses radix-2, mixed-radix, Rader, or Bluestein from transform structure;
- `KernelRadix2Plan::Auto` is an **experimental plan-time ISA selector** for one power-of-two radix-2 implementation family.

The latter measures scalar/AVX2/AVX-512 candidates during construction. Its tuning cost and selected ISA are part of the research record. It is not presented as a production-quality universal auto-tuner.

## Research frontier

The repository now covers one-dimensional, double-precision, single-threaded first-principles FFT algorithms, reusable complex/real plans, explicit x86 SIMD radix-2 codelets, experimental plan-time ISA selection, and a controlled serial FFTW baseline. Important next work includes:

- generated small-N and mixed-radix codelets;
- vectorized real-input recombination;
- cache-blocked, four-step, six-step, and cache-oblivious large transforms;
- Apple Accelerate/vDSP, Intel oneMKL, and pocketfft backends;
- Arm NEON/SVE codelets and physical x86-64/Arm baselines;
- hardware performance counters for cycles, instructions, cache misses, branches, and vector utilization;
- arbitrary-length reusable plans;
- `FFTW_PATIENT` / `EXHAUSTIVE` and wisdom persistence;
- multiple precisions and arbitrary-precision reference arithmetic;
- batched and multidimensional transforms;
- multithreaded scaling and GPU backends;
- persisted, statistically robust tuning wisdom rather than noisy cold per-plan selection.

## Research integrity

Headline results require raw data, exact source/binary/runtime provenance, sample/session counts, uncertainty analysis, and explicit timing semantics. SIMD claims additionally require runtime capability gates, explicit-width controls, documented persistent-state costs, and comparison against scalar plus each supported explicit ISA. Vendor comparisons require planner, normalization, allocation/alignment, representation, and threading semantics.

See [`docs/EXPERIMENTS.md`](docs/EXPERIMENTS.md), [`docs/SIMD_KERNELS.md`](docs/SIMD_KERNELS.md), [`docs/VENDOR_BENCHMARKS.md`](docs/VENDOR_BENCHMARKS.md), and [`CONTRIBUTING.md`](CONTRIBUTING.md).

## License

MIT. See [LICENSE](LICENSE).
