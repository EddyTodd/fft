# Development

The build is intentionally small. A top-level checkout builds the library, correctness tests, and example; when embedded with `add_subdirectory`, tests, examples, and install rules default off.

## Build and test

```bash
cmake --preset dev
cmake --build --preset dev
ctest --preset dev
```

Use the `release` preset for optimized builds.

Only three project options are needed:

- `FFTLAB_BUILD_TESTS` — build correctness tests;
- `FFTLAB_BUILD_EXAMPLES` — build the usage example;
- `FFTLAB_INSTALL` — enable install and `find_package` support.

Compiler warnings, sanitizers, coverage, and other developer policies are left to the parent project or ordinary CMake/compiler flags rather than wrapped in repository-specific helpers.

## Install

```bash
cmake --preset release
cmake --build --preset release
cmake --install build/release --prefix build/install
```

The installed package exports `fftlab::fftlab` for `find_package(fftlab CONFIG REQUIRED)`.

Performance measurement does not live in this repository. This project owns transform algorithms, plans, kernels, numerical correctness, and theory.
