# PR5 SIMD radix-2 kernel baseline

This evidence package studies machine-level optimization while holding the mathematical FFT family fixed. It compares the merged v3 reusable radix-2 plan, a new scalar codelet layout, explicit AVX2/FMA and AVX-512/FMA kernels, the experimental plan-time `Auto` selector, and FFTW `ESTIMATE` / `MEASURE` in the same randomized sessions.

## Formal matrix

- source commit: `31d901de63cd235ab13ad44a5f3210bdac0c0002`
- CPU model exposed by the container: AMD EPYC 9V74
- GCC 14.2, `-O3 -DNDEBUG -march=native -std=c++23`
- FFTW runtime: `fftw-3.3.10-sse2-avx`
- one thread
- power-of-two complex sizes: 64, 256, 1024, 4096, 16384, 65536
- 3 independently randomized sessions
- 31 execution samples per mode / size / session
- 1 setup sample per mode / size / session
- 7 policies
- **4,032 raw observations** total
- bootstrap analysis: 5,000 repetitions, seed `20260812`

Each execution observation times a persistent forward+inverse pair and divides elapsed time by two. Inverse normalization is included for every backend. Buffers and persistent plans are outside execution timing. Setup is measured separately; FFTW `MEASURE` setup forgets wisdom first.

## Modes

- `fftlab-plan / legacy` — merged v3 `Radix2Plan`
- `kernel / scalar` — compact swap list + stage-contiguous twiddles, scalar butterflies
- `kernel / avx2` — same plan structure with AVX2/FMA butterflies
- `kernel / avx512` — same plan structure with AVX-512/FMA butterflies
- `kernel / auto->ISA` — experimental plan-time empirical selection among supported kernels
- `fftw / estimate`
- `fftw / measure`

## Main findings

The best explicit SIMD codelet is **1.75–2.10x faster** than the merged v3 reusable plan over the formal size matrix and closes roughly **58–68%** of the latency gap from that plan to FFTW `MEASURE` in this environment.

Vector width has a real crossover. AVX2 is significantly faster at N=64, 256, and 1024; AVX-512 is significantly faster at N=4096, 16384, and 65536. Every AVX2-vs-AVX-512 bootstrap interval in the formal matrix excludes parity. This directly rejects a simplistic “always choose the widest supported ISA” policy.

The experimental `Auto` tuner is informative but not yet production-grade. It selects AVX2 consistently at small sizes and AVX-512 consistently at N=65536, but the selected ISA at N=1024–16384 varies across sessions and does not always match the pooled explicit-kernel winner. The complete 26–99 ms auto-tuning construction cost is reported rather than hidden.

Even after this optimization layer, the best local SIMD kernel remains approximately **1.97–2.35x slower than FFTW `MEASURE`** at the tested sizes. Therefore SIMD width alone does not explain FFTW's advantage; codelet generation, decomposition choice, scheduling, permutation strategy, cache behavior, and additional architecture-specific engineering remain research targets.

## Reproducibility

`metadata.json` records exact source Git blobs, compiler/build flags, formal binary SHA-256, raw-file SHA-256 and Git blob IDs, environment details, validation, and protocol semantics.

Regenerate the checked-in report directly from the compressed raw corpus:

```bash
python3 tools/analyze_kernel.py results/pr5-simd-kernel-baseline/raw \
  --bootstrap 5000 \
  --seed 20260812
```

The checked-in `ANALYSIS.md` was regenerated from those raw files and matched byte-for-byte during the final audit.

## Scope

These are descriptive results from one virtualized/containerized x86 environment. They do not establish a universal AVX2/AVX-512 crossover or a universal FFTW ranking. Physical-machine x86 and Arm studies, hardware counters, real-input SIMD post-processing, codelet generation, cache-blocked transforms, multithreading, and persisted tuning wisdom remain explicit future work.
