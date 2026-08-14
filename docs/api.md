# fftlab v1 API contract

## Precision and types

The stable v1 scalar concept accepts exactly `float` and `double`.

- `ComplexT<float>` / `VectorT<float>`: binary32 complex data.
- `ComplexT<double>` / `VectorT<double>`: binary64 complex data.
- `Complex32`, `Vector32`, `Complex64`, `Vector64`: convenience aliases.
- historical `Complex` and `Vector`: binary64 compatibility aliases.

## Mathematical convention

Forward:

`X[k] = sum_{n=0}^{N-1} x[n] exp(-2*pi*i*k*n/N)`.

Inverse uses the positive exponential and applies `1/N` normalization. All algorithm families use this same convention.

## Free algorithms

The free functions are primarily algorithm/reference interfaces. They may allocate temporary vectors and may construct roots/twiddles while executing.

Available mechanisms:

- `dft`
- `radix2`
- `recursive`
- `stockham`
- `radix4`
- `split_radix`
- `modified_split_radix`
- `mixed`
- `good_thomas`
- `rader`
- `bluestein`
- `transform`

`transform(input, Algo::Auto)` is a simple structural function dispatcher. Repeated workloads should prefer `Plan<T>`.

### Zero and one length

Free transforms accept `N=0` and return empty output. `N=1` is identity. Domain-specific algorithms still reject invalid nontrivial sizes (for example, radix-2 rejects non-power-of-two `N>1`).

## Reusable plans

### `Radix2Plan<T>`

Domain: power-of-two `N>=1`.

Persistent state: bit-reversal permutation and forward twiddles.

Execution: in-place or out-of-place; no caller scratch; no execution allocation/trigonometric setup.

### `RealRadix2Plan<T>`

Domain: power-of-two real `N>=1`.

Forward layout: `N/2+1` nonredundant complex bins (`1` bin for `N=1`).

Inverse consumes exactly that layout. Scratch requirement is `N/2` complex values for `N>1`.

### `SmallDftCodelet<T>`

Radices: 2, 3, 4, 5, 7.

Radix-2/4 use explicit butterflies. Radix-3/5/7 root matrices are precomputed at construction. Execution is allocation-free and trigonometry-free.

### `MixedRadixPlan<T>`

Reusable Cooley-Tukey plan using the small-codelet mechanism recursively. Top-level `Plan<T>` chooses it for supported smooth composite structures when PFA is not selected.

Scratch: `N` complex values. Execution allocates no dynamic memory and performs no trigonometric setup.

### `GoodThomasPlan<T>`

Domain: `N=a*b`, `gcd(a,b)=1`, `a,b>1`.

Persistent state includes CRT permutation maps plus reusable row/column plans. The Good-Thomas decomposition introduces no cross-stage twiddle multiplication (`twiddle_count()==0`).

Scratch: `N + 2*max(a,b)` complex values.

### `RaderPlan<T>`

Domain: prime `N>=3`.

Persistent state includes prime permutations, convolution spectrum, and reusable radix-2 convolution plan. If `N-1` is a power of two, the cyclic convolution is executed directly at that length.

### `BluesteinPlan<T>`

Domain: `N>=1`.

Persistent state includes chirps, convolution spectrum, and reusable radix-2 convolution plan.

### `Plan<T>`

Stable arbitrary-length structural planner. `PlanOptions` can force a supported mechanism or use `Structural` selection. The structural policy never times candidate algorithms.

`algorithm()` and `algorithm_name()` expose the actual selection. `plan_capabilities(N)` exposes candidate mechanism domains.

For out-of-place calls allocate `scratch_size()` complex elements. For in-place Rader/Bluestein calls allocate `inplace_scratch_size()` because those reductions require both algorithm scratch and temporary output.

`N=0/1` uses the identity plan.

## SIMD plan

`KernelRadix2Plan` is a binary64-only optional optimized extension. It accepts explicit `Scalar`, `Avx2`, or `Avx512` selection. Runtime capability checks prevent unsupported ISA execution. It performs no empirical tuning.

The portable f32/f64 planner is `Radix2Plan<T>`/`Plan<T>`.

## Aliasing

Use `*_inplace` APIs for explicit in-place execution. Out-of-place APIs are specified for distinct input/output storage. Scratch must not alias active input/output unless a specific function documents otherwise.

## Reuse and concurrency

Plans allocate/precompute only during construction. Execution is logically read-only with respect to plan state. A single plan object may be reused across threads **provided every concurrent call owns distinct mutable input/output/scratch storage**. v1 itself does not create threads or provide an internal thread pool.

## Correctness oracle

`oracle_dft` in `<fftlab/oracle.hpp>` computes the DFT in `long double` and is intended for deterministic subject-specific validation. `ErrorNorms`/`error_norms` provide normalized numerical error summaries. These are correctness utilities, not generic performance-statistics infrastructure.
