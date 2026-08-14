# Public API

## Headers

Use `<fftlab/fftlab.hpp>` for the complete installed API, or include focused headers directly.

- `types.hpp` — scalar concepts, complex/vector aliases, algorithm identifiers, and small number-theory helpers;
- `power2_algorithms.hpp` — direct DFT and power-of-two transform families;
- `arbitrary_algorithms.hpp` — mixed-radix, Good-Thomas, Rader, Bluestein, and generic dispatch;
- `plan.hpp`, `mixed_plan.hpp`, `good_thomas_plan.hpp`, `arbitrary_plan.hpp`, `convolution_plan.hpp`, `planner.hpp` — reusable structural plans;
- `codelet.hpp` — reusable small-radix codelets;
- `kernel.hpp` — explicit binary64 scalar/AVX radix-2 kernels;
- `oracle.hpp` — deterministic long-double DFT reference utilities.

## Scalar and container types

The library is precision-explicit:

```cpp
fftlab::ComplexT<float>
fftlab::ComplexT<double>
fftlab::VectorT<float>
fftlab::VectorT<double>

fftlab::Complex32
fftlab::Complex64
fftlab::Vector32
fftlab::Vector64
```

There are intentionally no unqualified historical `Complex` or `Vector` aliases in v1.

## Transform convention

`Direction::Forward` uses `exp(-2*pi*i*k*n/N)` and is unnormalized. `Direction::Inverse` uses the positive exponent and divides by `N`.

Most transform functions also provide a boolean convenience overload where `true` means inverse.

## Algorithm catalog and dispatch

`Algorithm` identifies the generic transform families. `all_algorithms` is the complete stable catalog.

- `algorithm_name(Algorithm)` returns the canonical CLI/report name;
- `parse_algorithm(name)` parses canonical names plus a small set of documented shorthand names;
- `supports_algorithm(algorithm, n)` reports structural domain support;
- `transform(input, algorithm, direction)` executes the selected family;
- `Algorithm::Auto` chooses a deterministic structural mechanism from the input length.

Unsupported explicit choices throw `std::invalid_argument`; they never silently fall back under the requested name.

The generic catalog currently includes direct DFT, radix-2 iterative/recursive, Stockham, radix-4, classical and modified split-radix, mixed-radix, Good-Thomas, Rader, and Bluestein.

## Reusable plans

Use plans for repeated transforms. Plans precompute structural choices, permutation data, twiddles/codelets, and scratch requirements so execution does not regenerate tables or allocate hidden work buffers.

`Plan<T>` is the deterministic arbitrary-length planner. It exposes the selected `PlanAlgorithm`, size, scratch requirement, and forward/inverse execution.

Specialized plans remain available when a caller wants an explicit mechanism:

- `Radix2Plan<T>`;
- `RealRadix2Plan<T>`;
- `MixedRadixPlan<T>`;
- `GoodThomasPlan<T>`;
- `RaderPlan<T>` / `BluesteinPlan<T>` through arbitrary/convolution planning;
- small-radix codelets.

See [`arbitrary-plans.md`](arbitrary-plans.md) and [`real-plans.md`](real-plans.md).

## Explicit SIMD kernels

`KernelRadix2Plan` is the binary64 fixed-width kernel plan. `KernelIsa` selects scalar, AVX2/FMA, or AVX-512/FMA. `kernel_capabilities()` reports runtime availability; requesting an unavailable ISA throws instead of silently substituting scalar execution.

The core library remains dependency-free and does not require FFTW or another vendor FFT package.

## Numerical reference

`oracle.hpp` provides a deterministic long-double DFT reference for correctness work. It is not a performance implementation.

## Research boundary

Algorithm crossover, setup-versus-reuse, planner policy quality, SIMD timing, vendor comparisons, replication, statistical inference, provenance, evidence, and reports belong in `EddyTodd/bench` and consume this public API.
