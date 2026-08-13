# fftlab

`fftlab` is a dependency-free C++23 library of one-dimensional CPU Fourier-transform algorithms and reusable plans. Version 1.0 freezes the repository around **algorithm mechanisms, structural planning, numerical correctness, and a reusable package API**. The benchmark campaigns accumulated during development remain in the repository only as legacy research material pending migration to [`EddyTodd/bench`](https://github.com/EddyTodd/bench).

## v1 at a glance

- complex binary32 (`std::complex<float>`) and binary64 (`std::complex<double>`);
- arbitrary transform lengths;
- sequential CPU execution;
- direct DFT, radix-2, recursive radix-2, Stockham, radix-4, classical split-radix, modified split-radix, mixed-radix, Good-Thomas/PFA, Rader, and Bluestein;
- reusable radix-2, mixed-radix, Good-Thomas, Rader, Bluestein, and structural arbitrary-length plans;
- real power-of-two plans with `N/2+1` Hermitian half spectra;
- reusable radix-2 small codelets for radices 2/3/4/5/7;
- explicit runtime-safe binary64 scalar/AVX2/AVX-512 radix-2 kernels;
- deterministic long-double DFT oracle utilities;
- installable/exported CMake package target `fftlab::fftlab`.

The precise scope and deferred domains are in [`docs/V1_SCOPE.md`](docs/V1_SCOPE.md).

## Transform convention

For `N` complex samples, the forward transform is

`X[k] = sum_n x[n] exp(-2*pi*i*k*n/N)`.

Forward transforms are unnormalized. Inverse transforms use the positive exponent and divide by `N`, so `inverse(forward(x))` returns `x` up to floating-point error.

The free algorithm functions treat `N=0` as an empty transform and `N=1` as identity. `Plan<T>` accepts both. Algorithm-specific plans that mathematically require a non-empty domain document and validate it explicitly.

## Precision architecture

The core API is templated only over the two supported v1 floating-point formats:

```cpp
fftlab::Vector32 x32;                   // vector<complex<float>>
fftlab::Vector64 x64;                   // vector<complex<double>>
fftlab::Plan<float>  p32(n);
fftlab::Plan<double> p64(n);
```

Historical `Complex` and `Vector` aliases remain binary64 compatibility aliases. The algorithm and planner implementations are shared templates rather than duplicated float/double source trees.

The x86 SIMD extension is intentionally binary64 in v1. Binary32 still has the full portable algorithm/planner catalog; dedicated f32 SIMD and Arm NEON are explicitly post-v1 work rather than hidden gaps in the portable API.

## Structural planner

`Plan<T>` is the stable arbitrary-length reusable planner. The default policy is deterministic and inspectable; it does **not** benchmark the host while constructing a plan.

Current structural policy:

1. `N=0/1` -> identity;
2. powers of two -> reusable radix-2;
3. composite lengths with a useful coprime split -> Good-Thomas/PFA by default;
4. remaining 2/3/5/7-decomposable composites -> planned mixed-radix;
5. rough composites -> Bluestein;
6. primes -> Rader only when its convolution is structurally shorter than Bluestein; otherwise Bluestein.

Every non-identity family can also be forced through `PlanOptions` when its domain is valid. `plan_capabilities(n)` exposes the structural possibilities for inspection/testing.

```cpp
#include <fftlab/fftlab.hpp>

fftlab::Plan<double> plan(60);
std::vector<std::complex<double>> input(60), output(60);
std::vector<std::complex<double>> scratch(plan.scratch_size());
plan.forward(input, output, scratch);
```

See [`docs/API.md`](docs/API.md) and [`docs/THEORY.md`](docs/THEORY.md).

## Planned mixed-radix and codelets

`MixedRadixPlan<T>` recursively plans Cooley-Tukey stages and uses reusable small DFT codelets for radices **2, 3, 4, 5, and 7**. Radix-2 and radix-4 have dedicated butterfly implementations; 3/5/7 root matrices are precomputed at codelet construction. Planned execution performs no dynamic allocation and no trigonometric setup.

For rough internal leaves the planner can precompute a direct DFT matrix, while the top-level structural `Plan<T>` normally chooses Bluestein instead of creating a large direct leaf.

## Good-Thomas / PFA

`GoodThomasPlan<T>` represents the coprime-factor prime-factor algorithm. It precomputes CRT input/output permutations, performs row/column reusable subplans, and requires **zero cross-stage Cooley-Tukey twiddle factors** at the Good-Thomas decomposition level. Explicit factor pairs are supported for auditability.

## Modified split-radix

`modified_split_radix()` implements a readable scaled conjugate-pair treatment following the Johnson-Frigo modified split-radix mechanism: recursive scale factors transform generic odd-branch twiddle multiplications into tangent/cotangent forms with fewer real multiplications. It is retained as an algorithmic reference mechanism rather than claiming generated-codelet-level arithmetic optimality.

Classical split-radix remains available separately.

## Real transforms

`RealRadix2Plan<T>` supports power-of-two real transforms in binary32 and binary64. Forward transforms return `N/2+1` nonredundant complex bins; inverse transforms consume that layout. Caller-owned complex scratch is explicit.

```cpp
fftlab::RealRadix2Plan<float> plan(1024);
std::vector<float> time(1024), restored(1024);
std::vector<std::complex<float>> spectrum(plan.spectrum_size());
std::vector<std::complex<float>> scratch(plan.scratch_size());
plan.forward(time, spectrum, scratch);
plan.inverse(spectrum, restored, scratch);
```

## SIMD extension

`KernelRadix2Plan` is an explicit binary64 extension with `Scalar`, `Avx2`, and `Avx512` modes. AVX modes are compiled behind function-level ISA targets on supported x86 compilers and are checked against runtime CPU capabilities before selection. Unsupported explicit requests throw; the scalar path always exists.

The old timing-based `Auto` kernel tuner is not part of the v1 core API. Empirical tuning belongs in the future benchmark/planner-research layer rather than being an implicit behavior of the only library path.

## Allocation, reuse, and thread safety

Plan construction may allocate persistent tables/permutations. After construction:

- `Radix2Plan<T>` and `KernelRadix2Plan` execute in-place without scratch allocation;
- `MixedRadixPlan<T>`, `GoodThomasPlan<T>`, `RaderPlan<T>`, `BluesteinPlan<T>`, and real plans use caller-owned scratch;
- planned execution does not allocate or regenerate trigonometric tables;
- plan execution methods are `const` and persistent plan state is not mutated.

A plan can therefore be reused concurrently when each invocation uses distinct input/output/scratch buffers. Use the explicit in-place APIs when aliasing is desired; out-of-place calls should use distinct input and output storage.

## Build, test, install

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
```

Sanitizers:

```bash
cmake -S . -B build-asan \
  -DCMAKE_BUILD_TYPE=Debug \
  -DFFTLAB_ENABLE_SANITIZERS=ON
cmake --build build-asan -j
ctest --test-dir build-asan --output-on-failure
```

Install and consume:

```bash
cmake --install build --prefix /your/prefix
```

```cmake
find_package(fftlab 1 CONFIG REQUIRED)
target_link_libraries(your_target PRIVATE fftlab::fftlab)
```

The package has no required FFTW dependency.

## Correctness strategy

The permanent deterministic suite covers both f32 and f64 and includes:

- an independent long-double `O(N^2)` DFT oracle;
- powers of two, mixed composites, PFA-friendly composites, primes, and awkward Bluestein fallbacks;
- zero/one/small lengths;
- forward/inverse round trips;
- impulse, constant, and single-tone identities;
- real-transform Hermitian symmetry;
- explicit codelet checks;
- scalar/AVX2/AVX-512 equivalence when those ISAs are available.

`<fftlab/oracle.hpp>` exposes the inexpensive long-double oracle for subject-specific correctness work. It is deliberately not a benchmark/statistics API.

## Repository layout

- `include/fftlab/` — stable installed API;
- `src/kernel.cpp` — compiled explicit x86 SIMD extension;
- `tests/` — deterministic v1 correctness tests;
- `docs/API.md` — detailed API/lifecycle contract;
- `docs/THEORY.md` — algorithm taxonomy/mechanisms;
- `docs/V1_SCOPE.md` — completeness boundary and deferred domains;
- `docs/LEGACY_RESEARCH.md` — development-era benchmark assets awaiting migration;
- `results/`, `tools/`, historical benchmark sources — legacy empirical research, not part of the installed library.

## Legacy research and `EddyTodd/bench`

The repository contains substantial historical timing campaigns, statistical analyzers, vendor comparisons, and raw result corpora from development. They remain useful evidence, but they are **not linked into `fftlab::fftlab` and are not installed**. Their intended destination is `EddyTodd/bench`.

See [`docs/LEGACY_RESEARCH.md`](docs/LEGACY_RESEARCH.md) for the concrete migration manifest.

## License

MIT. See [`LICENSE`](LICENSE).
