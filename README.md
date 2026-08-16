# fft

A C++23 library and study of **how different Fourier-transform algorithms actually work and compare** for sequential 1D CPU transforms.

The repository is useful independently: it contains reusable FFT/DFT algorithms, structural plans, SIMD kernels, deterministic numerical validation, and theory. Cross-algorithm performance research lives in [`EddyTodd/bench`](https://github.com/EddyTodd/bench).

## Use it

For a one-off transform:

```cpp
#include <fftlab/fftlab.hpp>

fftlab::Vector64 values{{1.0, 0.0}, {2.0, 0.0}, {3.0, 0.0}, {4.0, 0.0}};
const auto spectrum = fftlab::transform(values);
```

For repeated transforms, prefer `Plan<T>`. It uses `std::span` for input/output/scratch, so callers can use arrays, vectors, or other contiguous storage without adopting a library-specific container:

```cpp
fftlab::Plan<double> plan(1024);
std::vector<fftlab::Complex64> scratch(plan.scratch_size());
plan.forward(input, output, scratch);
```

As a CMake subdirectory:

```cmake
add_subdirectory(path/to/fft)
target_link_libraries(app PRIVATE fftlab::fftlab)
```

Or consume an installed package:

```cmake
find_package(fftlab 1 CONFIG REQUIRED)
target_link_libraries(app PRIVATE fftlab::fftlab)
```

## Build and test

```bash
cmake --preset release
cmake --build --preset release
ctest --preset release
```

A compile-checked plan example is included. The separate `package` preset validates installation, relocation, and a downstream `find_package` consumer without slowing the ordinary correctness loop.

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

`bench` studies direct-DFT/FFT crossovers, radix families, factorization effects, mixed-radix versus Good-Thomas, Rader versus Bluestein, setup versus reuse cost, planner quality, real transforms, SIMD benefit, and external baselines:

```bash
./bench run fft
```

## Read more

- [`docs/theory.md`](docs/theory.md) — mechanisms and derivations
- [`docs/api.md`](docs/api.md) — public API
- [`docs/arbitrary-plans.md`](docs/arbitrary-plans.md) — arbitrary-length planning
- [`docs/real-plans.md`](docs/real-plans.md) — real transforms
- [`docs/simd.md`](docs/simd.md) — architecture-specific kernels
- [`docs/development.md`](docs/development.md) — build options and package validation
- [`docs/scope.md`](docs/scope.md) — deliberate v1 boundary
- [`docs/references.md`](docs/references.md) — literature

## License

MIT.
