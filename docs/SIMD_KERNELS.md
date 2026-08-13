# SIMD and architecture-specific kernels

## v1 role

The stable portable FFT catalog is implemented for both binary32 and binary64. Architecture-specific codelets are an optional extension, not a requirement for correctness or planner availability.

`KernelRadix2Plan` provides explicit **binary64** x86 radix-2 execution modes:

- `Scalar`
- `Avx2` (requires AVX2 + FMA)
- `Avx512` (requires AVX-512F/DQ/VL + FMA)

The scalar path is always available.

## Runtime safety

On supported GCC/Clang x86 builds, AVX routines use function-level target attributes so compiling the library does not globally require the optional ISA. Runtime capability checks occur before an explicit AVX mode is accepted. Unsupported requests throw `std::invalid_argument`; unsupported instructions are never intentionally executed as feature probes.

## No implicit empirical tuner

The development-era kernel planner included a timing-based `Auto` mode. That mechanism is intentionally removed from the v1 core API. Host benchmarking, crossover inference, and persisted tuning wisdom belong to `EddyTodd/bench` or a future optional tuning layer.

The stable library path is deterministic: applications either use the portable `Plan<T>`/`Radix2Plan<T>` or explicitly request a supported binary64 kernel ISA.

## Binary32 and Arm

Binary32 is fully supported by the generic algorithms/plans, but v1 does not claim dedicated AVX binary32 kernels. Dedicated f32 SIMD is post-v1 optimization work.

Arm NEON is also explicitly deferred. The current completion environment cannot validate a NEON implementation on real Arm hardware, and an untested architecture-specific kernel would weaken rather than strengthen the portability contract. On Arm, the portable f32/f64 plans remain the v1 path. SVE is likewise post-v1.

## Numerical contract

FMA can change rounding because multiplication and addition may be fused. SIMD correctness is therefore tested against the generic radix-2 result using precision-appropriate tolerances rather than requiring bitwise identity. Forward and normalized inverse behavior must remain mathematically equivalent.

## Historical performance research

The previous SIMD timing corpus and auto-tuner study remain under `results/pr5-simd-kernel-baseline/` and related legacy sources/tools. Those assets are pending migration to `EddyTodd/bench`; they are not part of the installed v1 package.
