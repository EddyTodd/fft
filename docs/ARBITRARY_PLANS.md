# Reusable arbitrary-length FFT plans

## Research question

The legacy arbitrary-length APIs prove correctness and expose algorithm structure, but they rebuild chirps, permutations, convolution kernels, radix-2 twiddles, and temporary vectors on every call. That makes them inappropriate baselines for steady-state comparison with reused FFTW plans.

This milestone asks:

1. how much reusable precomputation changes Bluestein and Rader execution cost;
2. when planned Rader beats planned Bluestein for prime lengths;
3. how the structure of `p-1` changes Rader's convolution cost;
4. how quickly each plan repays setup relative to the legacy API;
5. how the first-principles planned reductions compare with persistent FFTW `ESTIMATE` and `MEASURE` plans under matched normalization/lifecycle semantics.

The central principle is the same as the earlier planning work: **setup cost and steady-state execution are different response variables**.

## Bluestein plan

For the forward DFT,

`X[k] = sum_n x[n] exp(-2 pi i n k / N)`.

Using

`2nk = n^2 + k^2 - (n-k)^2`,

Bluestein rewrites this as

`X[k] = c[k] * sum_n (x[n] c[n]) * conj(c[k-n]))`,

where

`c[j] = exp(-pi i j^2 / N)`.

`BluesteinPlan(N)` chooses

`M = next_pow2(2N - 1)`,

precomputes the N chirps, the M-point FFT of the symmetric convolution kernel, and an M-point `Radix2Plan`. Each execution then performs input chirp multiplication, one planned forward FFT, M pointwise products, one planned inverse FFT, and final chirp multiplication.

The caller supplies an M-complex scratch vector. Planned execution performs no dynamic allocation or trigonometric setup.

Inverse transforms use

`IDFT(x) = conj(DFT(conj(x))) / N`,

so one forward kernel spectrum is reused rather than storing a second kernel.

## Rader plan

For prime `p`, choose a primitive root `g` modulo `p`. Rader writes the nonzero DFT outputs as a cyclic convolution of length

`L = p - 1`:

`X[g^m] = x[0] + sum_q x[g^(-q)] exp(-2 pi i g^(m-q) / p)`.

`RaderPlan(p)` precomputes:

- output permutation `g^m mod p`;
- inverse input permutation `g^(-q) mod p`;
- forward convolution-kernel spectrum;
- reusable radix-2 convolution plan.

### Direct cyclic special case

When `L = p-1` is already a power of two, the cyclic convolution is evaluated directly with an L-point FFT. This differs materially from the legacy implementation, which always linearizes and zero-pads the convolution.

| p | Bluestein M | planned Rader M |
|---:|---:|---:|
| 17 | 64 | **16** |
| 257 | 1024 | **256** |
| 65537 | 262144 | **65536** |

When `p-1` is not a power of two, the current Rader implementation uses

`M = next_pow2(2(p-1) - 1)`

to evaluate a linear convolution and fold the `2L-1` result modulo L. A future reusable mixed-radix convolution planner may reduce this padding cost.

Inverse transforms use the same conjugation identity as Bluestein, so only one kernel spectrum is persistent.

## Evidence-derived `ArbitraryPlan` dispatch

The initial source freeze deliberately exposed both planned algorithms before changing automatic policy. The formal v6 matrix then showed:

- **N=17:** Rader 4.34× faster, with M=16 versus Bluestein M=64;
- **N=257:** Rader 4.33× faster, with M=256 versus Bluestein M=1024;
- when convolution lengths tie, differences are small: Rader wins N=31 by ~1.8%, Bluestein wins N=61, 127, 509, and 1009 by ~0.8–1.6%, and N=4093 is unresolved.

That evidence rejects a blanket `prime -> Rader` reusable policy. The final `ArbitraryPlan::Auto` rule is structural:

- power of two -> `Radix2Plan`;
- prime where `RaderPlan` has a **strictly shorter convolution** than `BluesteinPlan` -> `RaderPlan`;
- otherwise -> `BluesteinPlan`.

Explicit `Bluestein` and `Rader` policies remain available. The rule is intentionally conservative: it captures the large structural Rader advantage without baking a ~1% machine-specific tie-case crossover into the API. Future mixed-radix convolution planning or physical-hardware evidence may justify revising it.

## Buffer and memory contract

All complex values are `std::complex<double>` and caller buffers use ordinary `std::vector` storage.

### Bluestein persistent payload

- N chirp complex values;
- M kernel-spectrum complex values;
- M/2 radix-2 twiddles;
- M radix-2 permutation indices.

On a typical 64-bit ABI with 16-byte complex doubles and 8-byte `size_t`, the element payload is approximately

`16N + 32M bytes`,

excluding vector capacity/object/allocator overhead.

Execution requires one caller-owned M-complex scratch vector, approximately `16M` bytes.

### Rader persistent payload

- M kernel-spectrum complex values;
- M/2 radix-2 twiddles;
- 2(p-1) Rader permutation indices;
- M radix-2 permutation indices.

Approximate element payload:

`32M + 16(p-1) bytes`.

Execution requires one caller-owned M-complex scratch vector.

The extra persistent state is intentional: the experiment measures the memory-for-repeated-throughput tradeoff rather than treating precomputation as free.

## Correctness contract

The final specialized suite executes **5,574 arbitrary-plan checks**. It covers arbitrary, prime, composite, and power-of-two sizes; compares planned Bluestein and `ArbitraryPlan` with an independent direct DFT; exercises planned Rader on tested primes; verifies round trips; verifies domain rejection; and asserts the structural automatic policy.

When FFTW is available, `fft-arbitrary --self-test` additionally compares complete prime spectra at N=17, 127, and 509 for **653 FFTW frequency-bin cross-checks**.

The final policy source passed optimized GCC 14.2, optimized Clang 17, and GCC ASan/UBSan in the formal development environment.

## Formal benchmark semantics

Prime sizes:

`17, 31, 61, 127, 257, 509, 1009, 4093`.

Execution modes:

1. legacy Bluestein;
2. planned Bluestein;
3. legacy Rader;
4. planned Rader;
5. FFTW `ESTIMATE`;
6. FFTW `MEASURE`.

The legacy modes intentionally retain setup-inclusive behavior and are used to quantify the effect of planning. **Prime algorithm conclusions use planned-vs-planned execution.**

For persistent execution modes:

- plans and buffers exist before timing;
- forward and inverse are timed as a pair, then divided by two;
- inverse normalization is inside the timed path for every backend;
- caller input/output/scratch buffers are reused;
- execution-mode order is randomized inside each sample;
- size order is randomized independently per session;
- one common calibrated iteration count is used within each size;
- every raw observation is retained.

Setup samples are separate. FFTW `MEASURE` setup is cold: wisdom is forgotten immediately before each measured forward+inverse plan pair. FFTW arrays are allocated outside setup timing.

The formal corpus contains **4,560 observations**: 3 sessions × 8 sizes × (6 × 31 execution observations + 4 setup observations).

## Raw evidence transport

The canonical generated raw artifacts are gzip CSV streams and their SHA-256 hashes are recorded in `results/pr6-arbitrary-plan-baseline/metadata.json`.

This PR was written through a connector that can write UTF-8 files but not arbitrary binary contents, so each gzip stream is stored losslessly as two base64 text parts. `tools/analyze_arbitrary.py` detects `*.b64part*`, concatenates/decode/decompresses them in memory, and analyzes exactly the original CSV rows. The transport does not replace raw observations with aggregates.

## Predeclared hypotheses and outcome

**H8.** Persistent Bluestein and Rader plans materially outperform legacy setup-inclusive functions after a small number of reused transforms. **Supported:** Bluestein improves ~3.09–4.47×, Rader ~2.35–7.61×, with median setup break-even ~0.8–3.3 transforms.

**H9.** When `p-1` is a power of two, direct-cyclic planned Rader materially outperforms planned Bluestein because its convolution is shorter. **Supported** at N=17 and N=257, each by >4.3×.

**H10.** When both reductions use the same power-of-two convolution length, neither is assumed universally faster. **Supported:** differences are small and size-dependent, including one unresolved case.

**H11.** FFTW remains faster on most primes, but the gap narrows in structurally favorable Rader cases. **Supported:** FFTW `MEASURE` wins the matrix, while the best local plan is only ~1.28× slower at N=257 versus up to ~3.65× slower at N=4093.

## Threats to validity

- The formal host is virtualized; crossover claims are local to the recorded environment.
- Current non-direct Rader uses radix-2 zero-padded convolution when `p-1` is not a power of two. A mixed-radix convolution planner could materially change results.
- `std::vector<std::complex<double>>` is not an FFTW-aligned allocation contract.
- FFTW `MEASURE` explores plans beyond a fixed Rader/Bluestein decomposition, so it is a production reference rather than an algorithm-isolation control.
- Timing uses one deterministic random complex input family; numerical accuracy remains a separate axis.
- Setup amortization depends on allocator/cache state and the exact plan constructor.
- The structural automatic policy is intentionally conservative and should be re-evaluated after convolution-planner or hardware changes.

The result package records the exact formal source commit, post-evidence policy commit, binary hash, FFTW runtime, compiler/build, canonical gzip hashes, session protocol, and validation.
