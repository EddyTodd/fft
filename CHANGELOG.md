# Changelog

All notable changes to `fftlab` are recorded here. The installed C++ library follows semantic versioning; historical benchmark/research assets are not part of the permanent package surface.

## [Unreleased]

### Added

- Standard `dev`, `release`, `sanitize`, and checkout-local `package` CMake presets.
- External installed-package consumer smoke testing.
- Installed-package smoke now executes the linked scalar `KernelRadix2Plan` downstream after install/configure/build and verifies a deterministic forward/inverse transform.
- Build-time independent compilation of all 15 declared public FFT headers.
- Dependency-free release metadata verification across CMake, the public version header, citation metadata, and changelog history.

### Changed

- The installed header API is now an explicit target-owned CMake file set rather than directory-wide implicit installation.
- The MIT license is installed with the package.

## [1.0.0] - 2026-08-12

### Added

- Stable C++23 Fourier-transform library with direct DFT, iterative/recursive radix-2, Stockham, radix-4, classical split-radix, modified split-radix, mixed-radix, Good-Thomas, Rader, and Bluestein mechanisms.
- Reusable radix-2, real, mixed-radix, Good-Thomas, arbitrary-length Rader/Bluestein plans, and the general planner.
- Float and double support with scalar plus architecture-specific radix-2 kernel interfaces.
- Planned mixed-radix and arbitrary-length transforms with reusable scratch/workspace semantics.
- Deterministic correctness tests covering algorithm families, planners, kernels, round trips, and representative arbitrary sizes.
- Installable `fftlab::fftlab` CMake package and public version metadata.

### Architecture

- FFT algorithms, numerical semantics, planning, kernels, and correctness remain owned by `fft`.
- Generic timing/statistics/campaign/report/provenance machinery and cross-subject empirical analysis migrate to `EddyTodd/bench`.
