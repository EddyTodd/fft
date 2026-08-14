# fftlab

`fftlab` is a dependency-free C++23 library of sequential one-dimensional CPU Fourier transforms for binary32 and binary64. It provides representative transform mechanisms, reusable structural plans, SIMD kernels, deterministic numerical correctness, and an installable API. Empirical performance research lives in [`EddyTodd/bench`](https://github.com/EddyTodd/bench).

## Use

```cpp
#include <fftlab/fftlab.hpp>

std::vector<fftlab::Complex<double>> values = /* ... */;
fftlab::fft_radix2_iterative(values, fftlab::Direction::forward);
```

For repeated transforms, use a reusable plan so twiddles, structural choices, and scratch requirements are prepared once rather than regenerated per execution.

## Scope

Version 1 covers sequential 1D CPU transforms:

- direct DFT;
- iterative and recursive radix-2;
- Stockham radix-2;
- radix-4;
- classical and modified split-radix;
- mixed-radix Cooley-Tukey;
- Good-Thomas/PFA;
- Rader;
- Bluestein;
- reusable radix 2/3/4/5/7 codelets;
- radix-2, mixed-radix, Good-Thomas, Rader, Bluestein, real-radix2, and arbitrary plans;
- binary64 scalar/AVX2/AVX-512 radix-2 kernels;
- deterministic long-double DFT oracle utilities.

Multidimensional, batched, threaded, distributed, GPU, NUFFT, and empirical auto-tuning domains are outside v1.

## Convention

Forward transforms use `exp(-2*pi*i*k*n/N)` and are unnormalized. Inverse transforms use the positive exponent and divide by `N`.

## Build

```bash
cmake --preset dev
cmake --build --preset dev
ctest --preset dev
```

`release`, `sanitize`, and `package` presets use the same interface.

## Install

```bash
cmake --preset package
cmake --build --preset package
```

Consumer:

```cmake
find_package(fftlab 1 CONFIG REQUIRED)
target_link_libraries(your_target PRIVATE fftlab::fftlab)
```

The core library has no required FFTW or vendor dependency.

## Documentation

- [`docs/api.md`](docs/api.md) — public transform and plan API
- [`docs/theory.md`](docs/theory.md) — transform mechanisms and derivations
- [`docs/arbitrary-plans.md`](docs/arbitrary-plans.md) — arbitrary-length planning
- [`docs/real-plans.md`](docs/real-plans.md) — real-transform planning
- [`docs/simd.md`](docs/simd.md) — architecture-specific kernels
- [`docs/scope.md`](docs/scope.md) — v1 boundary
- [`docs/references.md`](docs/references.md) — literature

Crossover analysis, plan setup-versus-reuse studies, vendor comparisons, SIMD timing, statistics, provenance, and reports are intentionally centralized in [`EddyTodd/bench`](https://github.com/EddyTodd/bench).

## License

MIT.
