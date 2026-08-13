# Arbitrary-length reusable planning

## Stable API

The permanent arbitrary-length entry point is `Plan<T>` from `<fftlab/planner.hpp>`, with `T` equal to `float` or `double`.

`Plan<T>` separates a deterministic **structural policy** from any future empirical tuning. Construction records one explicit `PlanAlgorithm`, available through `algorithm()` / `algorithm_name()`.

## Structural policy

For the default `PlanOptions{}`:

1. `N=0/1` -> `Identity`;
2. power of two -> `Radix2`;
3. composite with a coprime factor split -> `GoodThomas` when `prefer_good_thomas=true`;
4. remaining 2/3/5/7-decomposable composite -> `MixedRadix`;
5. remaining composite -> `Bluestein`;
6. prime -> `Rader` only when its planned convolution is strictly shorter than Bluestein's, otherwise `Bluestein`.

`PlanPreference` can force radix-2, mixed-radix, Good-Thomas, Rader, or Bluestein when valid. Invalid forced choices throw.

## Planned mixed-radix

`MixedRadixPlan<T>` recursively decomposes `N=r*m`, with preferred small radices 4, 2, 3, 5, and 7. Stage twiddles and small-codelet roots are precomputed. Execution uses caller-owned `N`-complex scratch and performs no dynamic allocation or trigonometric setup.

## Good-Thomas / PFA

`GoodThomasPlan<T>` requires a coprime factorization `N=a*b`. The plan precomputes CRT input/output permutation maps, executes reusable row/column subplans, and introduces no cross-stage twiddle multiplication. `twiddle_count()` is therefore zero for the Good-Thomas decomposition itself.

## Rader

`RaderPlan<T>` requires prime `N>=3`. It precomputes primitive-root permutations and the FFT-domain convolution kernel. If `N-1` is a power of two, the cyclic convolution is executed directly at length `N-1`; otherwise a power-of-two zero-padded convolution is used and folded.

## Bluestein

`BluesteinPlan<T>` accepts any `N>=1`. It precomputes quadratic chirps, the FFT-domain convolution kernel, and a reusable radix-2 convolution plan.

## Scratch and in-place execution

`Plan<T>::scratch_size()` is the out-of-place requirement. Mixed/PFA/Rader/Bluestein expose their caller scratch explicitly.

Rader/Bluestein are naturally out-of-place reductions. `Plan<T>::inplace_scratch_size()` includes an additional `N` temporary output values for those algorithms; use that method when calling `forward_inplace` / `inverse_inplace`.

## Compatibility facade

`<fftlab/arbitrary_plan.hpp>` retains the historical binary64 `ArbitraryPlan` facade and explicit Rader/Bluestein policies. New code should use `Plan<T>` because it exposes the complete v1 planner and both precisions.

## Historical research

The prime crossover timing corpus and FFTW comparison under `results/pr6-arbitrary-plan-baseline/` are development evidence pending migration to `EddyTodd/bench`. They do not determine a hidden runtime tuner in v1.
