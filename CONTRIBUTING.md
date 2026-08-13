# Contributing to fftlab

fftlab's permanent role is a clear, reusable C++23 Fourier-transform algorithm library. Contributions should improve algorithms, planning, numerical correctness, portability, or package quality without re-coupling generic benchmark infrastructure to the core API.

## Algorithm contributions

A new permanent algorithm mechanism should have a clear reason to exist beyond increasing the catalog count. Document:

- transform convention and valid size domain;
- forward and inverse behavior;
- normalization;
- allocation/workspace requirements;
- asymptotic structure;
- primary literature reference for established algorithms;
- deterministic oracle and round-trip tests for f32/f64 where applicable.

## Plan contributions

Reusable plans must make lifecycle semantics explicit:

- persistent state built during construction;
- exact caller scratch requirement;
- in-place/out-of-place contract;
- whether execution allocates or evaluates trigonometric functions;
- plan reuse/concurrency expectations;
- structural policy/capability effects.

The default `Plan<T>` policy must remain deterministic and inspectable. Empirical timing/autotuning must be optional and separate rather than silently changing the only API path.

## Precision

The v1 scalar contract is binary32 and binary64. Shared templates are preferred over duplicated implementations. Precision-specific optimized code is acceptable when clearly isolated behind a portable fallback.

## SIMD

Architecture-specific kernels must:

- retain a portable scalar/generic path;
- runtime-gate optional ISAs before executing them;
- reject unsupported explicit requests;
- document exact ISA requirements;
- use numerical tolerance appropriate to FMA/rounding differences;
- never require optional ISA flags globally for the library.

Do not add an architecture-specific path that cannot be practically validated on the target architecture.

## Correctness

Run the deterministic suite, which covers the long-double DFT oracle, structural planner choices, f32/f64, real transforms, algorithm identities, and available explicit SIMD paths.

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DFFTLAB_WARNINGS_AS_ERRORS=ON
cmake --build build -j
ctest --test-dir build --output-on-failure
```

Sanitizers:

```bash
cmake -S . -B build-asan -DCMAKE_BUILD_TYPE=Debug \
  -DFFTLAB_ENABLE_SANITIZERS=ON -DFFTLAB_WARNINGS_AS_ERRORS=ON
cmake --build build-asan -j
ctest --test-dir build-asan --output-on-failure
```

Changes affecting packaging should also verify `cmake --install` and a separate `find_package(fftlab CONFIG REQUIRED)` consumer.

## Benchmark/research work

Do not add new generic statistics, campaign runners, vendor timing infrastructure, or raw campaign schemas to the installed core. The historical research layer is being consolidated into `EddyTodd/bench`; see [`docs/LEGACY_RESEARCH.md`](docs/LEGACY_RESEARCH.md).
