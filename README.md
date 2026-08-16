# fft

A C++23 study of **how different Fourier-transform algorithms actually work and compare** for sequential 1D CPU transforms.

This repository owns the FFT/DFT algorithms, reusable plans, numerical validation, and theory. Empirical performance research lives in [`EddyTodd/bench`](https://github.com/EddyTodd/bench).

## Try it

```cpp
#include <fftlab/fftlab.hpp>

fftlab::Vector64 values{{1.0, 0.0}, {2.0, 0.0}, {3.0, 0.0}, {4.0, 0.0}};
fftlab::radix2_inplace(values, fftlab::Direction::Forward);
```

Build and test:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

## Algorithms

The catalog represents important transform mechanisms rather than every FFT variant ever published:

- direct DFT;
- iterative and recursive radix-2;
- Stockham;
- radix-4;
- classical and modified split-radix;
- mixed-radix Cooley-Tukey;
- Good-Thomas / PFA;
- Rader;
- Bluestein;
- radix 2/3/4/5/7 codelets;
- reusable radix-2, mixed-radix, Good-Thomas, Rader, Bluestein, real, and arbitrary-length plans;
- binary64 scalar/AVX2/AVX-512 radix-2 kernels.

Planning stays here because decomposition, twiddles, scratch requirements, and structural choices are part of the FFT algorithm itself.

## Research

`bench` studies direct-DFT/FFT crossovers, radix families, factorization effects, mixed-radix versus Good-Thomas, Rader versus Bluestein, setup versus reuse cost, planner quality, real transforms, SIMD benefit, and external baselines such as FFTW when available.

From `bench`, the default study is intended to be:

```bash
python3 -m benchctl run fft
```

## Read more

- [`docs/theory.md`](docs/theory.md) — mechanisms and derivations
- [`docs/api.md`](docs/api.md) — public API
- [`docs/arbitrary-plans.md`](docs/arbitrary-plans.md) — arbitrary-length planning
- [`docs/real-plans.md`](docs/real-plans.md) — real transforms
- [`docs/simd.md`](docs/simd.md) — architecture-specific kernels
- [`docs/scope.md`](docs/scope.md) — deliberate v1 boundary
- [`docs/references.md`](docs/references.md) — literature

## License

MIT.
