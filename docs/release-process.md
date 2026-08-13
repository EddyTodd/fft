# Release process

An `fftlab` release freezes the reusable transform/planning API and its correctness semantics. Benchmark publication remains a separate `bench` concern.

## Version surfaces

Before tagging, synchronize:

- `project(fftlab VERSION ...)` in `CMakeLists.txt`;
- `include/fftlab/version.hpp`;
- `version` and `date-released` in `CITATION.cff`;
- the matching entry in `CHANGELOG.md`.

Run the dependency-free metadata verifier before the build matrix:

```sh
cmake -P cmake/VerifyReleaseMetadata.cmake
```

It derives the project version from CMake and fails if the public version header, citation version/date, or changelog release entry is inconsistent.

A release tag must point permanently at the exact validated source revision.

## Required validation

Run all standard standalone graphs from a clean checkout:

```sh
cmake --preset dev
cmake --build --preset dev
ctest --preset dev

cmake --preset release
cmake --build --preset release
ctest --preset release

cmake --preset sanitize
cmake --build --preset sanitize
ctest --preset sanitize

cmake --preset package
cmake --build --preset package
ctest --preset package
```

Every header in the `public_headers` file set must compile independently. The package graph must install `fftlab` and compile a separate `find_package(fftlab CONFIG REQUIRED)` consumer.

Changes to algorithms/plans must execute correctness coverage across representative power-of-two, composite, coprime-factor, prime, and awkward arbitrary lengths as applicable. Changes to SIMD/ISA kernels require real-target validation for every architecture claimed in release notes. Numerical claims must identify the norm/error contract and reference used.

## Planner/API compatibility

Before a release, review changes to:

- `Algo` names and selection semantics;
- plan construction and workspace requirements;
- forward/inverse normalization behavior;
- planner choice rules and capability reporting;
- reusable scratch/allocation guarantees;
- compatibility umbrella headers.

Behavioral changes that invalidate existing consumer assumptions require an appropriate semantic-version increment and explicit changelog entry.

## Research boundary

Do not delete vendor comparisons, kernel studies, arbitrary-length studies, real-plan studies, or historical empirical assets merely because generic measurement machinery exists in `bench`. Destructive migration requires treatment/correctness/evidence parity. Permanent FFT code must not depend on `bench`.

## Tag and GitHub release

After validation and merge:

1. create annotated `vMAJOR.MINOR.PATCH` at the exact validated `main` commit;
2. create a GitHub Release from that tag;
3. use the changelog entry as the release-note basis;
4. state architecture/vendor validation limitations explicitly;
5. verify rendered citation metadata matches the release version/date.

Never move a published release tag; publish a new version for corrections.
