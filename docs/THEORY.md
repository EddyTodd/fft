# FFT theory and algorithm taxonomy

## 1. Transform convention

For complex input \(x_0,\dots,x_{N-1}\), this project uses

\[
X_k = \sum_{n=0}^{N-1} x_n e^{-2\pi i kn/N}.
\]

The forward transform is unnormalized. The inverse uses the positive exponential and a final \(1/N\) scale. All implemented algorithms are intended to compute the same mathematical DFT; differences in measured output are floating-point effects, not different transform definitions.

## 2. Algorithms represented

| Implementation | Domain | Core idea | Time | Extra workspace | Research question |
|---|---|---|---:|---:|---|
| Direct DFT | any N | definition | \(O(N^2)\) | \(O(N)\) | baseline and reference cross-check |
| Iterative radix-2 DIT | powers of two | bit reversal + butterflies | \(O(N\log N)\) | \(O(1)\) besides returned copy | practical in-place baseline |
| Recursive radix-2 | powers of two | explicit divide-and-conquer | \(O(N\log N)\) | \(O(N)\) peak | recursion/allocation cost |
| Stockham radix-2 | powers of two | autosort/ping-pong stages | \(O(N\log N)\) | \(O(N)\) | regular access vs extra traffic |
| Radix-4 | powers of two | 4-way decomposition, radix-2 base | \(O(N\log N)\) | \(O(N)\) | fewer stages vs larger butterflies |
| Split-radix | powers of two | one N/2 + two N/4 transforms | \(O(N\log N)\) | \(O(N)\) | arithmetic count vs locality/overhead |
| Mixed-radix | composite N | Cooley–Tukey by factors | typically \(O(N\log N)\) | \(O(N)\) | factorization sensitivity |
| Rader | prime N | prime DFT → cyclic convolution | \(O(M\log M)\) | \(O(M)\) | prime-length alternative to Bluestein |
| Bluestein | any N | chirp-z → convolution | \(O(M\log M)\) | \(O(M)\) | robust arbitrary-length fallback |
| Auto | any N | inspectable dispatch policy | varies | varies | empirical dispatch research |

Here \(M\) is the power-of-two convolution size required by the reduction.

## 3. Structural operation models

`fft --complexity` reports a **structural complex-operation model**. It is intentionally not labeled as an exact CPU instruction count or FLOP count.

For a power-of-two radix-2 transform with \(L=\log_2 N\):

- complex butterfly additions/subtractions: \(N L\);
- nontrivial butterfly multiply slots: \((N/2)L\).

The implementation also performs twiddle-generation arithmetic, loop/index work, loads/stores, branches, allocation, and transcendental setup. Those costs are not represented by the structural model.

For split-radix, this repository uses the recurrence

\[
A(N)=A(N/2)+2A(N/4)+3N/2,
\]

\[
M(N)=M(N/2)+2M(N/4)+N/2,
\]

where multiplication by \(\pm i\) is treated as trivial. This model is useful for comparing decomposition structure, not for claiming an exact real-arithmetic record.

The 2007 modified split-radix work of Johnson and Frigo reduces the real arithmetic count below classical split-radix. More recent work has lowered theoretical leading constants further. This repository therefore treats “fewest operations” and “fastest implementation” as separate research dimensions.

## 4. Prime-length transforms

### Rader

For prime \(N\), nonzero indices form a cyclic multiplicative group. A primitive root permutes those indices, transforming the DFT into a cyclic convolution of length \(N-1\). This implementation computes that convolution using zero-padding plus radix-2 FFTs.

### Bluestein

Bluestein rewrites the phase using a quadratic identity and turns an arbitrary-size DFT into a convolution. It is particularly useful when Cooley–Tukey factorization is unfavorable.

Rader and Bluestein have similar asymptotic forms but different convolution lengths, setup work, data permutations, and numerical behavior. Their crossover should be measured, not assumed.

## 5. Why operation count is insufficient

Modern performance also depends on:

- cache and TLB behavior;
- memory traffic and access regularity;
- branch and instruction throughput;
- SIMD/vectorization opportunities;
- allocation strategy;
- twiddle generation/storage;
- compiler transformations;
- transform size and factorization;
- in-place vs out-of-place semantics;
- setup/planning amortization.

This is why the project reports theoretical models beside empirical timings rather than collapsing both into a single ranking.

## 6. Numerical error

Floating-point FFTs are not exact. Important mechanisms include rounding in butterflies, inaccurate twiddles, error accumulation through recursive/staged decomposition, and cancellation. The benchmark therefore reports normalized forward and backward errors in \(L_1\), \(L_2\), and \(L_\infty\), plus a conventional forward→inverse round-trip maximum error.

The current reference is a long-double \(O(N^2)\) DFT, not arbitrary-precision arithmetic. This is much stronger than comparing algorithms only with each other, but it is still a documented limitation. A future milestone should add MPFR or another independently controlled high-precision oracle as an optional research dependency.
