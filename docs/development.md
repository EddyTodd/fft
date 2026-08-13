# Development workflow

`fftlab` uses the shared subject-repository developer workflow.

## Standard presets

```sh
cmake --preset dev
cmake --build --preset dev
ctest --preset dev
```

Optimized configuration:

```sh
cmake --preset release
cmake --build --preset release
ctest --preset release
```

Strict local sanitizer configuration:

```sh
cmake --preset sanitize
cmake --build --preset sanitize
ctest --preset sanitize
```

The sanitizer preset enables `FFTLAB_ENABLE_SANITIZERS=ON` and `FFTLAB_WARNINGS_AS_ERRORS=ON`.

## Research boundary

The shared presets build and test the permanent v1 library only. Historical research C++ sources under `research/` and Python campaign/analysis programs under `tools/` are migration assets, not part of the normal development graph.

New empirical orchestration belongs in `EddyTodd/bench`; FFT-specific algorithm, planner, numerical-correctness, and kernel work remains here.

## User-local configuration

Put machine/compiler/SDK overrides in `CMakeUserPresets.json`; it is intentionally ignored by Git. Shared presets must remain portable and must not encode a local vendor FFT installation or ISA-specific workstation assumptions.
