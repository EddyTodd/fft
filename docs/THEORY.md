# FFT theory and v1 algorithm taxonomy

## Transform convention

For complex input `x[0..N-1]`, fftlab computes

`X[k] = sum_n x[n] exp(-2*pi*i*k*n/N)`.

Forward transforms are unnormalized. Inverse transforms use the positive exponential and apply `1/N` normalization.

## Representative catalog

| Mechanism | Domain | Key role in v1 |
|---|---|---|
| Direct DFT | any N | mathematical definition and small-size oracle comparison |
| Iterative radix-2 | powers of two | canonical in-place Cooley-Tukey baseline |
| Recursive radix-2 | powers of two | explicit divide-and-conquer form |
| Stockham radix-2 | powers of two | autosort/ping-pong staging |
| Radix-4 | powers of two | larger radix butterfly family |
| Classical split-radix | powers of two | N/2 + two N/4 decomposition |
| Modified split-radix | powers of two | Johnson-Frigo scaled conjugate-pair arithmetic reduction mechanism |
| Mixed-radix | composite | Cooley-Tukey over nonuniform factors |
| Good-Thomas / PFA | coprime composite factors | CRT permutation with no cross-stage twiddles |
| Rader | prime N | prime DFT -> cyclic convolution |
| Bluestein | any N | chirp-z reduction -> convolution |
| Small codelets | radix 2/3/4/5/7 | planned small-factor kernels |

The free algorithm functions emphasize readable mechanisms. Reusable plans precompute state and separate construction from execution.

## Planned mixed-radix

For `N=r*m`, a Cooley-Tukey decomposition recursively transforms `m`-point subsequences, applies stage twiddles, and combines `r` values with a small DFT. `MixedRadixPlan<T>` stores the stage twiddles and reusable codelet state during construction.

The v1 codelet set covers radices 2, 3, 4, 5, and 7. Radix-2/4 use explicit butterflies. Radix-3/5/7 use precomputed root matrices, so plan execution does not evaluate trigonometric functions.

## Good-Thomas / Prime Factor Algorithm

When `N=a*b` with `gcd(a,b)=1`, the Chinese Remainder Theorem can map the one-dimensional DFT index into a two-dimensional product group. The transform then becomes `a`- and `b`-point transforms plus input/output permutations, without the cross twiddle factors required by ordinary Cooley-Tukey at that decomposition level.

`GoodThomasPlan<T>` explicitly stores the CRT permutations and exposes `twiddle_count()==0` for the top-level PFA decomposition. Its row/column transforms are themselves reusable mixed-radix plans.

## Classical and modified split-radix

Classical split-radix decomposes a power-of-two DFT into one even `N/2` transform and two odd `N/4` transforms. It is retained directly because it is historically and structurally important.

The modified split-radix family reduces real arithmetic further by recursively rescaling subtransforms. In the Johnson-Frigo construction, these scale factors turn selected ordinary complex twiddle products into tangent/cotangent forms that require fewer real multiplications. `modified_split_radix()` implements this scaled conjugate-pair mechanism in a readable recursive form.

This implementation is a faithful mechanism/reference treatment, not a claim that a readable recursive C++ routine attains the exact instruction count of generated codelets or is the fastest machine implementation.

## Rader and Bluestein

Rader maps a prime-length DFT to cyclic convolution of length `N-1`. The reusable plan directly uses an `N-1` radix-2 cyclic convolution when `N-1` is already a power of two; otherwise it uses a zero-padded linear convolution and folds the result.

Bluestein maps any DFT length to convolution using quadratic chirps. Its reusable plan stores chirps and the FFT-domain convolution kernel.

The structural planner compares convolution lengths for primes: Rader is selected only when it produces a strictly shorter planned convolution; otherwise Bluestein is the conservative default.

## Real FFT reduction

For even power-of-two `N=2M`, a real sequence can be packed into one `M`-point complex sequence by placing even samples in the real part and odd samples in the imaginary part. A post-processing identity reconstructs the `N/2+1` unique frequency bins implied by Hermitian symmetry. `RealRadix2Plan<T>` implements this reduction for f32/f64.

## Planning is part of the algorithm interface

A reusable FFT library must distinguish mathematical decomposition from reusable state. v1 plans expose scratch sizes and perform no dynamic allocation/trigonometric table generation during execution. Structural dispatch is deterministic and inspectable; empirical timing/autotuning is deliberately outside the stable planner.

## Numerical behavior

Floating-point FFTs accumulate rounding error through butterflies, twiddle multiplication, cancellation, and recursive/staged composition. Binary32 and binary64 therefore have separate numerical contracts. The deterministic suite compares both against a long-double direct DFT oracle and validates round-trip/identity properties.

Long double is a practical deterministic oracle, not arbitrary precision. Optional higher-precision research can be added later without changing the v1 transform API.
