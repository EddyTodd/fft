# fftlab

`fftlab` is a dependency-free C++23 library of one-dimensional CPU Fourier-transform algorithms and reusable plans for binary32 and binary64.

This repository owns transform mechanisms, structural planning, numerical correctness, SIMD implementation, and theory. Empirical timing, vendor comparisons, campaign orchestration, statistics, provenance, and reports live in [`EddyTodd/bench`](https://github.com/EddyTodd/bench).

## v1 catalog

- direct DFT;
- iterative and recursive radix-2;
- Stockham radix-2;
- radix-4;
- classical and modified split-radix;
- mixed-radix Cooley-Tukey;
- Good-Thomas / PFA;
- Rader;
- Bluestein;
- reusable small codelets for radices 2/3/4/5/7;
- reusable radix-2, mixed-radix, Good-Thomas, Rader, Bluestein, real-radix2, and arbitrary structural plans;
- explicit binary64 scalar/AVX2/AVX-512 radix-2 kernels;
- deterministic long-double DFT oracle utilities.

The arbitrary-length `Plan<T>` chooses a deterministic structural mechanism and exposes capabilities/options for explicit family selection. Plan execution uses caller-owned scratch where needed and avoids per-execution table regeneration/allocation.

See `docs/API.md`, `docs/THEORY.md`, and `docs/V1_SCOPE.md`.

## Transform convention

Forward transforms use `exp(-2*pi*i*k*n/N)` and are unnormalized. Inverse transforms use the positive exponent and divide by `N`.

## Build and test

```bash
cmake --preset dev
cmake --build --preset dev
ctest --preset dev
```

Release, sanitizer, and package presets are also provided. The permanent deterministic suite covers f32/f64, arbitrary structural domains, forward/inverse identities, real transforms, codelets, and scalar/AVX equivalence where supported.

## Install

```bash
cmake --preset package
cmake --build --preset package
```

```cmake
find_package(fftlab 1 CONFIG REQUIRED)
target_link_libraries(your_target PRIVATE fftlab::fftlab)
```

The core has no required FFTW dependency.

## Scope

Version 1 is sequential 1D CPU Fourier transforms. Multidimensional/batched/threaded/distributed/GPU transforms, DCT/DST, NUFFT, sparse FFTs, generated codelet systems, additional SIMD architectures, and empirical auto-tuning are post-v1 domains.

Algorithm crossover, plan setup-versus-reuse, real-transform performance, explicit SIMD timing, and FFTW/vendor comparisons are research questions for `EddyTodd/bench`, using this stable library API.

## License

MIT.
