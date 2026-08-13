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

## Install-package check

The `package` preset installs a Release build into `build/package-prefix`:

```sh
cmake --preset package
cmake --build --preset package
ctest --preset package
```

The install exports `fftlab::fftlab`, all permanent public FFT/planning headers, relocatable package/version files, and the project license. The 15 supported installed headers are an explicit `public_headers` CMake file set, so an internal or experimental header cannot become installed API merely by appearing under `include/fftlab/`. The prefix is local to the checkout and can be consumed by a separate project through `CMAKE_PREFIX_PATH` or `fftlab_DIR`.

Normal non-sanitized standalone CTest graphs also run `fftlab.package-consumer`. The test installs the current build into an isolated prefix, configures a separate project with `find_package(fftlab 1 CONFIG REQUIRED)`, and compiles it against the installed `fftlab::fftlab` target. The consumer includes all 15 installed headers, so the manifest, exported target, and installed tree must agree. Sanitizer configurations omit this downstream distribution smoke because sanitizer runtime requirements are a separate consumer contract.

## Research boundary

The shared presets build and test the permanent v1 library only. Historical research C++ sources under `research/` and Python campaign/analysis programs under `tools/` are migration assets, not part of the normal development graph.

New empirical orchestration belongs in `EddyTodd/bench`; FFT-specific algorithm, planner, numerical-correctness, and kernel work remains here.

## User-local configuration

Put machine/compiler/SDK overrides in `CMakeUserPresets.json`; it is intentionally ignored by Git. Shared presets must remain portable and must not encode a local vendor FFT installation or ISA-specific workstation assumptions.
