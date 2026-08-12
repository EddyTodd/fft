# SIMD radix-2 kernel research

This document defines the v5 research layer for studying how much of the gap between a readable radix-2 FFT and a production FFT library can be explained by **data layout, reusable codelet state, SIMD width, fused arithmetic, and lightweight empirical dispatch**.

The goal is not to replace the algorithm taxonomy with ISA-specific code. The goal is to hold the mathematical decomposition fixed while changing one machine-level implementation dimension at a time.

## 1. Research questions

The kernel experiment asks:

1. How much performance is gained by replacing a full bit-reversal table with a precomputed list of actual swap pairs?
2. Does storing twiddles contiguously per FFT stage improve the scalar implementation enough to justify additional persistent plan memory?
3. How much additional speed comes from two-complex-wide AVX2/FMA butterflies?
4. How much additional speed comes from four-complex-wide AVX-512/FMA butterflies?
5. Is AVX-512 consistently faster than AVX2 merely because the vectors are twice as wide?
6. Can a small plan-time empirical tuner choose the faster supported kernel without embedding a machine-specific static table?
7. How much of the previously measured FFTW execution gap is closed by this one optimization layer, and what gap remains?
8. How many transforms are required to repay the larger codelet plan or the auto-tuning step?

These are implementation hypotheses. Neither ISA support nor vector width is treated as proof of superior performance.

## 2. Mathematical transform remains radix-2 Cooley-Tukey

`KernelRadix2Plan` computes the same decimation-in-time power-of-two radix-2 transform as the existing `Radix2Plan`.

For each butterfly,

```text
u = x[k]
v = w * x[k + N/2]
x[k]       = u + v
x[k + N/2] = u - v
```

where `w` is the stage twiddle. Forward transforms use negative phase; inverse transforms conjugate the stored forward twiddle and apply `1/N` normalization.

Keeping the decomposition fixed is important: scalar-vs-AVX2-vs-AVX-512 comparisons are intended to study code generation and vector width, not different FFT algorithms.

## 3. Swap-list permutation

The v3 reusable plan stores one bit-reversed index for every input position and tests `i < bit_reverse[i]` during every execution.

The kernel plan instead stores only pairs that actually need to be exchanged:

```text
(i, reverse_bits(i)) where i < reverse_bits(i)
```

Execution therefore walks a compact swap list without a per-index branch.

This is still an explicit bit-reversal permutation. The milestone does not claim it is the best possible permutation strategy; Stockham/autosort, out-of-place codelets, and cache-blocked permutations remain separate research directions.

## 4. Stage-contiguous twiddles

The v3 plan retains `N/2` global twiddles and accesses them with a stage-dependent stride.

The kernel plan stores each stage's twiddle sequence contiguously. Across all radix-2 stages this requires

```text
1 + 2 + 4 + ... + N/2 = N - 1
```

complex twiddles.

The extra storage removes strided twiddle lookup inside the butterfly loop and lets a SIMD load fetch consecutive twiddles directly.

This is an explicit memory-for-throughput tradeoff. It must be charged to plan storage and setup rather than described as a free optimization.

## 5. SIMD representation

`Complex` remains `std::complex<double>`; the public buffer type does not change. The x86 kernels use unaligned vector loads/stores so callers are not required to provide a new aligned allocator.

On supported GCC/Clang x86 builds:

- AVX2 processes two complex doubles per 256-bit vector;
- AVX-512 processes four complex doubles per 512-bit vector;
- both use FMA for the real/imaginary complex multiply recombination;
- the inverse normalization loop is vectorized at the same width.

Small stage tails remain scalar when there are fewer butterflies than one full vector.

The functions are compiled with function-level target attributes, so the surrounding binary does **not** need to require AVX2 or AVX-512 globally. Runtime CPUID support is checked before a kernel can be selected.

On compilers/architectures where these function-level x86 targets are not compiled, the scalar kernel remains available and explicit unsupported SIMD requests are rejected.

## 6. FMA and numerical semantics

A fused multiply-add performs a multiply and add/subtract with one final rounding instead of the intermediate rounding of two separate operations. Therefore an FMA kernel is not expected to be bit-identical to the scalar implementation even when it is mathematically equivalent.

Correctness tests use a tight norm-scaled tolerance against the existing radix-2 implementation and independently verify forward/inverse round trips. The benchmark never treats bitwise identity as a prerequisite for a valid floating-point FFT implementation.

The broader numerical-accuracy framework remains in `docs/EXPERIMENTS.md`; this milestone focuses on machine-level performance while preserving correctness.

## 7. Runtime capabilities

`kernel_capabilities()` reports whether the binary/CPU combination can safely execute the explicit research kernels.

The current x86 requirements are:

- AVX2 path: AVX2 + FMA;
- AVX-512 path: AVX512F + AVX512DQ + AVX512VL + FMA.

Capability only means **legal to execute**. It does not mean that the wider path is faster.

## 8. Lightweight empirical auto-tuning

`KernelIsa::Auto` deliberately does not choose the widest supported ISA.

After plan state is constructed, the tuner:

1. creates a deterministic complex input;
2. evaluates every supported candidate: scalar, AVX2, AVX-512;
3. runs five timing rounds;
4. rotates/reverses candidate order across rounds to reduce systematic order bias;
5. uses a size-scaled iteration count;
6. selects the candidate with the lowest median forward/inverse execution time;
7. records the total tuning time in `tuning_ns()`.

The auto tuner is intentionally much lighter than FFTW `MEASURE`, but it is **not free**. Formal results report its complete plan-construction cost and the selected ISA per session.

The explicit AVX2 and AVX-512 modes remain available so the tuner can be audited against the same raw data.

## 9. Persistent memory model

Ignoring container/allocator overhead and stage-offset bytes:

### v3 `Radix2Plan`

- `N` `size_t` bit-reversal entries;
- `N/2` complex twiddles.

On a conventional 64-bit ABI this is approximately `16N` bytes of retained element payload.

### `KernelRadix2Plan`

- one pair of `size_t` values per actual bit-reversal swap;
- `N-1` complex stage-local twiddles;
- `log2(N)` stage offsets.

The exact swap count depends on N. The twiddle payload alone is approximately `16(N-1)` bytes, so the SIMD-friendly plan is intentionally larger than the v3 plan.

A speedup that depends on larger reusable state must therefore be interpreted together with setup latency and persistent memory.

## 10. Benchmark matrix

`fft-kernel` compares seven persistent complex-transform policies under one harness:

1. `fftlab-plan / legacy` — v3 reusable `Radix2Plan`;
2. `kernel / scalar` — swap-list + stage-contiguous twiddles, scalar butterflies;
3. `kernel / avx2` — same plan structure, AVX2/FMA butterflies;
4. `kernel / avx512` — same plan structure, AVX-512/FMA butterflies;
5. `kernel / auto->ISA` — one empirically selected persistent kernel;
6. `fftw / estimate` — persistent FFTW `ESTIMATE` plans;
7. `fftw / measure` — persistent FFTW `MEASURE` plans.

Execution semantics match the v4 vendor study:

- one thread;
- double precision;
- power-of-two complex transforms;
- persistent plans;
- buffers outside the timed region;
- forward + inverse pairs with elapsed time divided by two;
- inverse normalization included for every backend;
- randomized mode order inside every sample;
- randomized transform-size order across formal sessions;
- setup measured separately;
- FFTW `MEASURE` setup measured cold after forgetting wisdom;
- every raw sample retained.

This matrix allows the v3-to-SIMD improvement and the remaining SIMD-to-FFTW gap to be measured in the **same sessions**.

## 11. Statistical interpretation

`tools/analyze_kernel.py` reports:

- pooled median execution latency;
- independent bootstrap 95% speedup intervals;
- AVX2-vs-AVX-512 uncertainty;
- auto-selection frequency by session;
- best explicit SIMD speedup over the v3 plan;
- fraction of the v3-to-FFTW-`MEASURE` latency gap closed;
- plan setup medians;
- setup break-even transform counts.

A vector-width winner is called only when the bootstrap interval for the pairwise speedup excludes parity.

## 12. Known limitations

This milestone is intentionally narrow:

- x86 SIMD codelets only; Arm NEON/SVE is not yet implemented;
- power-of-two complex transforms only;
- real-input post-processing is not yet vectorized;
- bit reversal itself remains scalar;
- no cache blocking or six-step/four-step decomposition;
- no generated small-N codelet family;
- no explicit prefetching study;
- no hardware-counter data yet;
- auto-tuning has no persisted wisdom/cache across processes;
- compiler target-attribute support currently provides SIMD kernels on GCC/Clang x86; other compiler backends remain future work;
- the formal baseline is still a virtualized machine, not a controlled physical-hardware study.

These limitations define the next experiments rather than being hidden implementation debt.
