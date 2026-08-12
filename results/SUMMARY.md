# Research results summary

## v5 — SIMD radix-2 codelet baseline

The v5 milestone asks how much of the gap between a readable reusable radix-2 plan and FFTW can be explained by **plan-state layout and explicit SIMD codelets while the mathematical decomposition stays fixed**.

The checked-in `pr5-simd-kernel-baseline/` corpus contains **4,032 raw observations** from 3 randomized sessions at N=64, 256, 1024, 4096, 16384, and 65536. The exact source commit, source blobs, formal binary hash, compiler/runtime, raw hashes, and benchmark semantics are recorded in metadata.

### Execution

| N | v3 plan | Scalar codelet | AVX2/FMA | AVX-512/FMA | FFTW MEASURE | Best SIMD / v3 |
|---:|---:|---:|---:|---:|---:|---:|
| 64 | 313.6 ns | 294.4 ns | **179.4 ns** | 215.0 ns | 82.1 ns | **1.747×** |
| 256 | 1.60 µs | 1.53 µs | **758.1 ns** | 936.4 ns | 367.5 ns | **2.104×** |
| 1024 | 7.62 µs | 7.30 µs | **3.75 µs** | 4.41 µs | 1.70 µs | **2.035×** |
| 4096 | 41.13 µs | 38.22 µs | 24.46 µs | **21.07 µs** | 10.69 µs | **1.952×** |
| 16384 | 228.53 µs | 202.27 µs | 142.49 µs | **121.85 µs** | 51.76 µs | **1.876×** |
| 65536 | 1.227 ms | 993.61 µs | 733.51 µs | **642.54 µs** | 292.70 µs | **1.909×** |

The best explicit SIMD path is **1.747–2.104× faster** than the merged v3 plan and closes roughly **57.9–68.2%** of the latency gap from that plan to FFTW `MEASURE`. FFTW still remains approximately **1.97–2.35× faster** than the best local SIMD codelet across the matrix.

### Wider vectors are not universally faster

The recorded data have a clear size-dependent crossover:

- AVX2 is faster at N=64, 256, and 1024;
- AVX-512 is faster at N=4096, 16384, and 65536;
- all six AVX2-versus-AVX-512 bootstrap 95% intervals exclude parity.

This is an implementation/environment result rather than a universal x86 rule. It demonstrates why ISA availability or vector width alone is not a defensible dispatch policy.

### Layout matters before vectorization

The scalar codelet uses the same mathematical radix-2 decomposition but changes reusable plan state:

- only actual bit-reversal swap pairs are stored;
- each stage's twiddles are stored contiguously;
- the execution loop therefore removes a per-index permutation branch and strided twiddle lookup.

The scalar path is already faster than the v3 plan at the tested sizes, though the major gain comes from explicit SIMD/FMA.

### Planning economics

The SIMD-friendly plan stores more twiddle state than the v3 plan. Its explicit setup cost remains small enough to amortize quickly in this baseline:

- about **20 transforms** at N=64;
- about **9.6** at N=256;
- about **2.4** at N=1024;
- about **1.35** at N=4096;
- below **one transform** at N=16384 and 65536 according to the recorded median setup/execution deltas.

The experimental `Auto` selector is different. It times scalar/AVX2/AVX-512 candidates during construction and costs roughly **26–99 ms** in this experiment. Its choices are not perfectly stable: it consistently selects AVX2 at small sizes and AVX-512 at N=65536, but at intermediate sizes its session-level choice can vary and can disagree with the pooled explicit-kernel winner.

That mismatch is preserved as evidence, not hidden. `Auto` is a research mechanism for studying adaptive planning; it is not presented as a production-quality universal dispatch policy.

### Remaining FFTW gap

One vectorized radix-2 codelet closes most of the original local-to-FFTW latency difference, but not all of it. The remaining ~2× gap points toward additional engineering dimensions including:

- generated/specialized codelet families;
- decomposition and schedule selection;
- permutation strategy;
- cache blocking and large-transform algorithms;
- alignment/data-layout choices;
- real-input vectorization;
- persisted planning wisdom;
- other architecture-specific optimization.

See `pr5-simd-kernel-baseline/ANALYSIS.md` and `docs/SIMD_KERNELS.md`.

## v4 — FFTW production-library baseline

The v4 corpus contains **3,456 raw observations** and established that persistent FFTW execution is much faster than the v3 fftlab plans while planner economics can reverse the end-to-end winner for short workloads.

Across its formal matrix, FFTW `MEASURE` was approximately **3.77–4.68× faster** than fftlab for planned complex execution and **4.16–6.59× faster** for planned real execution. Cold `MEASURE` planning cost tens of milliseconds to roughly two seconds and required approximately **19,668 to 2,126,344 transforms** to amortize versus `ESTIMATE`, depending on size and workload.

The methodological conclusion was that steady-state execution speed and total workload cost are separate questions.

## v3 — reusable planning and real-input specialization

The v3 corpus contains **2,790 raw observations**. It established that the same radix-2 decomposition becomes **1.295–1.593× faster** when reusable permutation/twiddle setup is removed from steady-state execution, with plan construction amortizing after roughly **1.75–3.95 transforms**. The specialized real path reaches roughly **1.74–1.75×** the planned complex throughput at the larger tested sizes.

## v2 — algorithm families, structural complexity, and accuracy

Selected results from the formal expanded-algorithm study include:

- Rader about **1.29× faster than Bluestein** at N=509 and N=1009 in the recorded environment;
- iterative radix-2 about **1089× faster than direct DFT** at N=1024;
- a pedagogical split-radix implementation about **4.51× slower** than iterative radix-2 at N=256 despite a lower structural multiplication model;
- substantial algorithm-dependent forward-error differences across five deterministic signal families.

That milestone demonstrated that mathematical operation count, numerical behavior, and machine latency are related but not interchangeable.

## Evidence scope

All current formal baseline numbers come from virtualized/containerized development environments. They support implementation- and methodology-specific conclusions, not universal hardware rankings. Physical x86-64 and Arm runs with controlled power/affinity/thermal state are required before making architecture-general crossover claims.

The v1 baseline remains under `baseline-linux-amd-epyc-gcc14.csv` for historical continuity.
