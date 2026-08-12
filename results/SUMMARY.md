# Research results summary

## v6 — reusable arbitrary-length plans and prime crossover

The v6 milestone asks whether the repository's earlier prime-FFT ranking survives once **Bluestein and Rader both receive reusable planning and matched lifecycle semantics**.

The checked-in `pr6-arbitrary-plan-baseline/` evidence represents **4,560 raw observations** from 3 randomized sessions at N=17, 31, 61, 127, 257, 509, 1009, and 4093. The exact formal source commit, source blobs, binary hash, compiler/runtime, original gzip hashes, and benchmark semantics are recorded in metadata.

### Planning is a first-order effect

| N | Legacy Bluestein | Planned Bluestein | Legacy Rader | Planned Rader |
|---:|---:|---:|---:|---:|
| 17 | 2.66 µs | 744.9 ns | 1.31 µs | 171.6 ns |
| 31 | 3.43 µs | 767.1 ns | 2.55 µs | 753.2 ns |
| 61 | 7.10 µs | 1.65 µs | 5.26 µs | 1.66 µs |
| 127 | 14.96 µs | 3.63 µs | 11.04 µs | 3.68 µs |
| 257 | 52.43 µs | 16.96 µs | 23.72 µs | 3.92 µs |
| 509 | 66.75 µs | 17.39 µs | 49.48 µs | 17.67 µs |
| 1009 | 142.14 µs | 38.02 µs | 110.21 µs | 38.62 µs |
| 4093 | 680.10 µs | 216.11 µs | 512.46 µs | 218.05 µs |

Planning improves Bluestein by roughly **3.09–4.47×** and Rader by roughly **2.35–7.61×** relative to their historical setup-inclusive APIs. Median plan construction generally repays itself within about **0.8–3.3 transforms**.

### Rader's real structural advantage is convolution length

For prime `p`, Rader reduces the nonzero DFT to a cyclic convolution of length `p-1`. The planned implementation evaluates that convolution directly when `p-1` is already a power of two.

| N | Bluestein convolution M | Rader convolution M | Planned winner | Bluestein / Rader median ratio |
|---:|---:|---:|---|---:|
| 17 | 64 | **16 cyclic** | **Rader** | 4.340× |
| 31 | 64 | 64 | **Rader** | 1.018× |
| 61 | 128 | 128 | **Bluestein** | 0.992× |
| 127 | 256 | 256 | **Bluestein** | 0.986× |
| 257 | 1024 | **256 cyclic** | **Rader** | 4.327× |
| 509 | 1024 | 1024 | **Bluestein** | 0.984× |
| 1009 | 2048 | 2048 | **Bluestein** | 0.984× |
| 4093 | 8192 | 8192 | unresolved | 0.991× |

At N=17 and N=257, Rader is more than **4.3× faster** because it uses a convolution one quarter the Bluestein length. When the power-of-two convolution lengths tie, the differences are small: Rader wins N=31 by ~1.8%; Bluestein wins N=61, 127, 509, and 1009 by roughly 0.8–1.6%; and N=4093 has a bootstrap interval spanning parity.

This supports a structural `ArbitraryPlan::Auto` rule: **use Rader for a prime only when it produces a strictly shorter convolution; otherwise use Bluestein**. Explicit policies remain available so this rule can be challenged by future hardware/mixed-radix evidence.

### Why this revises—not contradicts—the v2 result

The v2 study found legacy Rader about **1.29× faster than legacy Bluestein** at N=509 and N=1009. Those APIs allocate and rebuild reusable state on every call. The v6 study measures persistent plans and finds planned Bluestein slightly faster at both sizes.

Both results are valid for their stated API semantics. The changed ranking demonstrates the repository's central methodology point: **planning/setup semantics can change an apparent algorithm winner**.

### Production baseline

FFTW `MEASURE` remains faster at every formal prime size. The best local planned reduction is about **1.28× slower at N=257** and up to about **3.65× slower at N=4093** in this environment. The narrow N=257 gap is itself informative: a structurally favorable direct-cyclic Rader reduction can approach an adaptive production library much more closely than the generic zero-padded reductions.

See `pr6-arbitrary-plan-baseline/ANALYSIS.md` and `docs/ARBITRARY_PLANS.md`.

## v5 — SIMD radix-2 codelet baseline

The v5 milestone asks how much of the gap between a readable reusable radix-2 plan and FFTW can be explained by **plan-state layout and explicit SIMD codelets while the mathematical decomposition stays fixed**.

The checked-in `pr5-simd-kernel-baseline/` corpus contains **4,032 raw observations** from 3 randomized sessions at N=64, 256, 1024, 4096, 16384, and 65536.

| N | v3 plan | Scalar codelet | AVX2/FMA | AVX-512/FMA | FFTW MEASURE | Best SIMD / v3 |
|---:|---:|---:|---:|---:|---:|---:|
| 64 | 313.6 ns | 294.4 ns | **179.4 ns** | 215.0 ns | 82.1 ns | **1.747×** |
| 256 | 1.60 µs | 1.53 µs | **758.1 ns** | 936.4 ns | 367.5 ns | **2.104×** |
| 1024 | 7.62 µs | 7.30 µs | **3.75 µs** | 4.41 µs | 1.70 µs | **2.035×** |
| 4096 | 41.13 µs | 38.22 µs | 24.46 µs | **21.07 µs** | 10.69 µs | **1.952×** |
| 16384 | 228.53 µs | 202.27 µs | 142.49 µs | **121.85 µs** | 51.76 µs | **1.876×** |
| 65536 | 1.227 ms | 993.61 µs | 733.51 µs | **642.54 µs** | 292.70 µs | **1.909×** |

The best explicit SIMD path is **1.747–2.104× faster** than the merged v3 plan and closes roughly **57.9–68.2%** of its latency gap to FFTW `MEASURE`. FFTW remains approximately **1.97–2.35× faster** than the best local SIMD codelet.

AVX2 is faster at N=64, 256, and 1024; AVX-512 is faster at N=4096, 16384, and 65536 in the formal host. All six pairwise bootstrap intervals exclude parity. The experimental `Auto` tuner is deliberately retained as imperfect: its 26–99 ms construction cost and occasional disagreement with the pooled explicit winner are research results, not hidden failures.

## v4 — FFTW production-library baseline

The v4 corpus contains **3,456 raw observations** and established that persistent FFTW execution is much faster than the v3 fftlab plans while planner economics can reverse the end-to-end winner for short workloads.

Across its formal matrix, FFTW `MEASURE` was approximately **3.77–4.68× faster** than fftlab for planned complex execution and **4.16–6.59× faster** for planned real execution. Cold `MEASURE` planning cost tens of milliseconds to roughly two seconds and required approximately **19,668 to 2,126,344 transforms** to amortize versus `ESTIMATE`, depending on size and workload.

## v3 — reusable planning and real-input specialization

The v3 corpus contains **2,790 raw observations**. It established that the same radix-2 decomposition becomes **1.295–1.593× faster** when reusable permutation/twiddle setup is removed from steady-state execution, with plan construction amortizing after roughly **1.75–3.95 transforms**. The specialized real path reaches roughly **1.74–1.75×** the planned complex throughput at the larger tested sizes.

## v2 — algorithm families, structural complexity, and accuracy

Selected results from the setup-inclusive expanded-algorithm study include:

- legacy Rader about **1.29× faster than legacy Bluestein** at N=509 and N=1009;
- iterative radix-2 about **1089× faster than direct DFT** at N=1024;
- a pedagogical split-radix implementation about **4.51× slower** than iterative radix-2 at N=256 despite a lower structural multiplication model;
- substantial algorithm-dependent forward-error differences across five deterministic signal families.

That milestone demonstrated that mathematical operation count, numerical behavior, and machine latency are related but not interchangeable. v6 further demonstrates that lifecycle semantics are another independent axis.

## Evidence scope

All current formal baseline numbers come from virtualized/containerized development environments. They support implementation- and methodology-specific conclusions, not universal hardware rankings. Physical x86-64 and Arm runs with controlled power/affinity/thermal state are required before architecture-general claims.

The v1 baseline remains under `baseline-linux-amd-epyc-gcc14.csv` for historical continuity.
