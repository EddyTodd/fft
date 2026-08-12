# Planning, execution semantics, and real-input FFTs

This document defines the v3 research model for separating **one-time transform setup** from **steady-state execution** and for exploiting the conjugate symmetry of real-input DFTs.

The distinction is necessary before comparing this repository with FFTW, Accelerate/vDSP, oneMKL, pocketfft, or other production libraries. A benchmark that charges one implementation for planning and twiddle construction while another reuses precomputed state is not measuring equivalent work.

## 1. Why a plan exists

The original iterative radix-2 implementation performs two kinds of work on every transform:

1. structural FFT work: permutation and butterfly evaluation;
2. reusable setup work: stage-root evaluation and generation of the twiddle sequence.

`Radix2Plan` moves reusable information out of the execution path. Construction precomputes:

- the complete bit-reversal permutation for the selected transform size;
- `N/2` forward twiddle factors `exp(-2πik/N)`.

Inverse execution conjugates the stored forward twiddles and applies the usual `1/N` normalization. The execution methods neither resize buffers nor allocate workspace and do not call `sin`, `cos`, or another trigonometric function.

This gives the repository two intentionally different radix-2 measurements:

- **legacy complex** — the existing in-place radix-2 kernel, including repeated stage-root/twiddle construction;
- **planned complex** — the same radix-2 decomposition with permutation/twiddles precomputed once.

The comparison isolates the practical value of reusable planning without changing the mathematical FFT family.

## 2. Plan cost and amortization

Planning is not free. For `Radix2Plan(N)` the retained state is:

- `N` permutation indices;
- `N/2` complex twiddles.

The current constructor computes each bit-reversed index explicitly and evaluates each stored twiddle independently. Its construction cost is therefore an implementation property, not a theoretical lower bound on FFT planning.

The benchmark records plan construction separately. A descriptive break-even count is

```text
complex plan break-even = median(plan setup) /
                          [median(legacy execution) - median(planned execution)]
```

when the denominator is positive.

This answers a practical question: **how many transforms of the same size must be executed before this plan has repaid its construction cost?**

The result is machine/compiler/implementation specific. It must not be treated as a universal FFT constant.

## 3. Real-input transform representation

For real input `x[n]`, the DFT satisfies Hermitian symmetry:

```text
X[N-k] = conj(X[k]).
```

Only `N/2 + 1` complex bins are unique for even `N`: DC (`k=0`), positive frequencies, and the Nyquist bin (`k=N/2`). `RealRadix2Plan` therefore exposes a packed half-spectrum of exactly that size rather than returning an N-element complex array containing redundant conjugates.

The current specialization supports power-of-two `N >= 2`.

## 4. N/2-complex reduction

Let `N = 2M`. Pack the even and odd real samples into one M-point complex sequence:

```text
z[j] = x[2j] + i x[2j+1],    0 <= j < M.
```

Compute the M-point complex FFT `Z[k]`. Define

```text
A[k] = 1/2 ( Z[k] + conj(Z[M-k]) )
B[k] = -i/2 ( Z[k] - conj(Z[M-k]) ).
```

`A[k]` is the M-point transform of the even samples and `B[k]` is the M-point transform of the odd samples. The N-point real-input DFT is then

```text
X[k] = A[k] + exp(-2πik/N) B[k],    0 <= k <= M.
```

Indices are interpreted modulo `M` where required. This produces the `M+1 = N/2+1` unique bins directly.

The inverse reconstructs the packed complex spectrum from the half-spectrum, performs one M-point inverse complex FFT, then unpacks real and imaginary parts back into even and odd time samples.

## 5. Why the real path can be faster

A full N-point complex FFT applied to real-valued data ignores known structure. The specialized path instead performs:

- one N/2-point planned complex FFT;
- O(N) packing/recombination work;
- half-spectrum storage rather than a redundant N-bin complex output.

The asymptotic class remains `O(N log N)`, but the constant factors and memory traffic are reduced. The exact speedup is empirical because packing, post-processing, cache behavior, compiler optimization, and vectorization all matter.

Using the repository's canonical radix-2 structural model with `M=N/2`, the complex core contributes approximately

```text
M log2(M) complex additions
(M/2) log2(M) complex multiplications
```

plus O(N) recombination. The current forward recombination performs roughly three complex add/subtract operations and one nontrivial twiddle multiplication per returned bin before accounting for endpoint/trivial-factor simplifications. This is intentionally a structural description rather than a CPU-FLOP claim.

The specialization should therefore be evaluated on two axes:

1. **structural advantage** — fewer nonredundant frequency values and a smaller complex FFT;
2. **measured advantage** — observed latency and memory behavior on a concrete machine.

## 6. Persistent state and caller memory

Planning trades repeated setup arithmetic for retained metadata. The storage model is explicit so a latency win cannot hide an unreported memory cost.

### Complex radix-2 plan

`Radix2Plan(N)` retains:

```text
N     size_t indices
N/2   complex<double> twiddles
```

On a conventional 64-bit target where `sizeof(size_t)=8` and `sizeof(complex<double>)=16`, the element payload is approximately

```text
8N + 16(N/2) = 16N bytes
```

excluding vector/object headers, allocator bookkeeping, and spare capacity.

Planned complex execution is in-place and requires no plan-owned scratch buffer. The caller's transform array itself is `N` complex doubles, or approximately `16N` bytes.

### Real radix-2 plan

`RealRadix2Plan(N)` contains an `M=N/2` complex plan plus `M+1` post-processing twiddles. Its retained element payload is therefore approximately

```text
M size_t indices                 = 4N bytes
(M/2) complex twiddles           = 4N bytes
(M+1) complex post-twiddles      = 8N + 16 bytes
------------------------------------------------
total                            ≈ 16N + 16 bytes
```

again excluding object/container/allocator overhead.

For one real forward transform the caller supplies:

```text
N doubles input                  = 8N bytes
(N/2+1) complex half-spectrum    = 8N + 16 bytes
N/2 complex scratch              = 8N bytes
------------------------------------------------
concurrent caller payload        ≈ 24N + 16 bytes
```

The inverse has the same leading-order caller payload: half-spectrum input, N real outputs, and N/2 complex scratch.

A full N-bin complex spectrum would require `16N` bytes by itself; the packed real spectrum requires approximately `8N+16`, asymptotically halving frequency-domain output storage. The scratch buffer means total working-set behavior still requires empirical measurement rather than inferring cache performance from output size alone.

These byte estimates are representation-level models. `sizeof` values and allocator overhead should be recorded directly when publishing results on unusual ABIs or alternative complex representations.

## 7. Execution benchmark contract

`fft-plan` uses forward+inverse pairs so the same preallocated buffer can be reused indefinitely without copying a fresh input into the timed region. Pair elapsed time is divided by two and reported as nanoseconds per transform.

For each transform size the benchmark records five modes:

- `complex-setup` — construct a fresh `Radix2Plan`;
- `real-setup` — construct a fresh `RealRadix2Plan`;
- `legacy-complex` — forward+inverse pairs using the existing in-place radix-2 kernel;
- `planned-complex` — forward+inverse pairs using one reused `Radix2Plan`;
- `planned-real` — forward+inverse pairs using one reused `RealRadix2Plan` and preallocated half-spectrum/scratch buffers.

Important semantics:

- planned execution contains no allocation and no trigonometric setup;
- legacy and planned complex paths both operate in place;
- real output and scratch storage are allocated before timing;
- calibration uses the planned-complex path and a common iteration count is then used for the execution modes;
- mode order is randomized independently inside every sample;
- formal experiments additionally randomize transform-size order across sessions;
- all raw samples are retained.

## 8. Planned-real break-even

A second break-even metric asks whether the real specialization's setup overhead is justified relative to a reusable N-point complex plan:

```text
real break-even = max(0, median(real setup) - median(complex setup)) /
                  [median(planned complex execution) - median(planned real execution)]
```

This is intentionally conservative: if the real plan costs no more to construct than the full complex plan in a recorded run, the extra setup penalty is zero.

## 9. Correctness requirements

The dedicated planned self-test checks:

- planned complex forward transforms against the existing radix-2 implementation;
- planned complex forward/inverse round trips;
- the packed real half-spectrum against the nonredundant bins of an N-point complex FFT of the same real input;
- real forward/inverse round trips;
- domain rejection for unsupported transform sizes.

The planned API does not resize user buffers. Incorrect buffer sizes are rejected rather than silently allocating or truncating data.

## 10. What this enables next

The plan abstraction establishes the benchmark contract required for fair production-library comparisons. Future backends should expose, or be normalized to, at least these phases:

1. plan/setup construction;
2. required persistent plan memory;
3. required caller/workspace memory;
4. execution with already-created plans and already-allocated buffers;
5. destruction, excluded from execution latency.

Production-library studies must additionally control:

- in-place vs out-of-place transforms;
- real vs complex input;
- forward normalization convention;
- planning effort/flags;
- thread count;
- precision;
- alignment requirements;
- transform size and batch count;
- cold-plan latency vs repeated steady-state throughput.

Without these controls, a headline library ranking would be scientifically weak even if the timing code itself were precise.

## 11. Current limitations

The v3 plan layer is deliberately narrow:

- plans currently specialize power-of-two radix-2 transforms;
- real-input plans currently require power-of-two even lengths;
- no SIMD-specialized codelets are generated;
- plan construction itself is not optimized aggressively;
- no wisdom/cache file is persisted across processes;
- no multithreaded planning or execution is implemented;
- the benchmark does not yet collect hardware performance counters.

These constraints make the experiment interpretable. Subsequent milestones can relax them one at a time while preserving the setup/execution separation introduced here.
