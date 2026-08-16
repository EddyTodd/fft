# Development

A top-level checkout builds the library, deterministic numerical tests, and a small plan example. When `fftlab` is embedded with `add_subdirectory`, tests, examples, and install rules default off so the dependency does not alter its parent's build surface.

## Normal workflow

```bash
cmake --preset dev
cmake --build --preset dev
ctest --preset dev
```

Use `release` for optimized builds and `sanitize` for ASan/UBSan with warnings-as-errors.

Important options:

- `FFTLAB_BUILD_TESTS` — deterministic FFT/planner correctness tests;
- `FFTLAB_BUILD_EXAMPLES` — small standalone usage examples;
- `FFTLAB_INSTALL` — install/export/package rules (top-level by default);
- `FFTLAB_ENABLE_SANITIZERS` — ASan/UBSan on supported non-MSVC compilers;
- `FFTLAB_WARNINGS_AS_ERRORS` — promote compiler warnings to errors.

## Installed package validation

Package validation is intentionally separate from the normal test loop:

```bash
cmake --preset package
cmake --build --preset package
ctest --preset package
```

The package preset enables `FFTLAB_BUILD_PACKAGE_TESTS`; its package-consumer test installs into an isolated prefix, relocates it, and validates a separate `find_package(fftlab CONFIG REQUIRED)` consumer. Ordinary `dev` and `release` tests do not pay that cost.

Machine/compiler/SDK overrides belong in ignored `CMakeUserPresets.json`; shared presets should remain portable and should not encode workstation-specific vendor FFT or ISA paths.

Performance experiments belong in `EddyTodd/bench`; transform algorithms, planning, numerical correctness, and kernels remain here.
