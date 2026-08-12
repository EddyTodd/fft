# v4 FFTW production-library baseline

This dataset is the repository's first controlled comparison against an optimized external FFT implementation. It compares the reusable fftlab plans introduced in v3 with the system FFTW 3.3.10 runtime under matched complex/real and setup/execution semantics.

## Protocol

- source commit: `781ad3deca0d49afaa9010fe01077f0dd0d11fa4`
- FFTW runtime: `fftw-3.3.10-sse2-avx`
- compiler: GCC 14.2, `-O3 -DNDEBUG -march=native -std=c++23`
- one thread
- sizes: 64, 256, 1024, 4096, 16384, 65536
- 3 independently randomized sessions
- 31 execution samples per backend/planner/kind/session
- 1 independently cold setup sample per backend/planner/kind/session
- 5 warmups
- 2 ms execution calibration target

Every execution datum is the elapsed time of a forward+inverse pair divided by two. FFTW's inverse is unnormalized, so the required `1/N` scaling is included inside its timed execution path. All transform buffers are allocated before execution timing. Setup timing excludes caller-buffer allocation.

FFTW `MEASURE` setup is deliberately cold: FFTW wisdom is forgotten before every measured forward+inverse plan pair. `ESTIMATE` and `MEASURE` execution both use persistent plans.

## Raw evidence

`raw/` contains three gzip-compressed CSV session shards. The corpus contains exactly **3,456 observations**. `metadata.json` records SHA-256 and Git blob hashes for each raw file plus the exact source blobs and formal binary hash.

Regenerate the report with:

```bash
python3 tools/analyze_vendor.py results/pr4-fftw-baseline/raw \
  --seed 20260812 \
  --output results/pr4-fftw-baseline/ANALYSIS.md
```

## Interpretation

The dataset answers two different questions:

1. **steady-state execution:** how fast is an already-planned transform?
2. **amortized workload cost:** after including plan creation, which strategy is cheapest for K repeated transforms of the same size?

The distinction matters. FFTW `MEASURE` is the fastest steady-state path throughout this dataset, but its cold planning cost can require tens of thousands to millions of transforms to amortize relative to FFTW `ESTIMATE`. At smaller sizes fftlab's much cheaper plan can minimize total cost for short workloads despite slower execution.

These numbers describe this FFTW binary and a virtualized/containerized AMD EPYC environment. They are not universal library rankings. Future external-validity work should repeat the protocol on named physical x86-64 and Arm machines, multiple FFTW builds, and platform libraries such as Accelerate and oneMKL.
