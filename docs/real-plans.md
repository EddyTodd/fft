# Reusable power-of-two and real FFT plans

## `Radix2Plan<T>`

`Radix2Plan<float>` and `Radix2Plan<double>` are the basic reusable power-of-two plans. Construction precomputes bit-reversal indices and forward twiddles. Forward/inverse execution is in-place and requires no caller scratch.

Inverse execution conjugates the stored forward twiddles and applies `1/N` scaling.

## `RealRadix2Plan<T>`

For power-of-two `N>=2`, the real plan packs even/odd samples into one `N/2`-point complex transform and uses the standard recombination identity to return only the `N/2+1` nonredundant frequency bins. `N=1` is supported as an identity special case.

Forward:

- input: `N` real values;
- output: `N/2+1` complex values (`1` for `N=1`);
- scratch: `N/2` complex values for `N>1`.

Inverse consumes that exact half-spectrum layout and returns `N` real values.

## Lifecycle contract

Plan construction may allocate and evaluate trigonometric functions. Execution does not allocate or regenerate twiddle tables. Persistent plan state is immutable during execution; concurrent reuse is safe when each call owns distinct mutable buffers.

## Numerical validation

The deterministic v1 suite validates both f32 and f64 real plans against the long-double direct DFT oracle, checks every stored half-spectrum bin against the corresponding full complex transform, verifies Hermitian symmetry, and performs forward/inverse round trips.

## Historical benchmark material

The earlier plan/reuse timing experiment is retained under `results/pr3-planning-real-baseline/` and related tools, but it is legacy empirical research pending migration to `EddyTodd/bench`. It is not part of the installed library contract.
