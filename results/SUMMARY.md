# Research results summary

## v4 — FFTW production-library baseline

The v4 milestone asks a different question from earlier algorithm races: **which FFT strategy is cheapest for a workload once planning and execution are both modeled?**

The checked-in `pr4-fftw-baseline/` corpus contains **3,456 raw observations** from 3 randomized sessions at N=64, 256, 1024, 4096, 16384, and 65536. It compares fftlab reusable complex/real plans with FFTW 3.3.10 `ESTIMATE` and `MEASURE`. Source/runtime/build hashes and exact benchmark semantics are recorded in metadata.

### Steady-state execution

| N | fftlab complex | FFTW ESTIMATE | FFTW MEASURE | fftlab real | FFTW ESTIMATE real | FFTW MEASURE real |
|---:|---:|---:|---:|---:|---:|---:|
| 64 | 313.7 ns | 137.3 ns | 81.6 ns | 305.5 ns | 85.6 ns | 66.8 ns |
| 256 | 1.58 µs | 468.3 ns | 370.8 ns | 1.35 µs | 323.6 ns | 204.5 ns |
| 1024 | 7.77 µs | 2.13 µs | 1.66 µs | 6.03 µs | 1.03 µs | 918.6 ns |
| 4096 | 40.55 µs | 15.01 µs | 10.77 µs | 27.35 µs | 5.13 µs | 4.72 µs |
| 16384 | 230.84 µs | 108.18 µs | 51.97 µs | 138.47 µs | 41.48 µs | 25.11 µs |
| 65536 | 1.20 ms | 354.77 µs | 287.35 µs | 711.37 µs | 211.74 µs | 171.15 µs |

Across this matrix, FFTW `MEASURE` is approximately **3.77–4.68× faster** than fftlab for planned complex execution and **4.16–6.59× faster** for planned real execution. The relevant bootstrap 95% intervals remain well above parity.

### Planning changes the answer

Cold planning is not free. In the recorded run:

- FFTW `ESTIMATE` complex setup ranges from roughly 92 µs to 738 µs;
- FFTW `MEASURE` complex setup ranges from about 49 ms to 1.87 s;
- FFTW `MEASURE` real setup ranges from about 21 ms to 2.05 s.

As a result, the faster `MEASURE` execution path may require **19,668 to 2,126,344 repeated transforms** to repay its additional planning cost relative to `ESTIMATE`, depending on transform size and real/complex workload.

The same effect appears when comparing fftlab with FFTW `ESTIMATE` at small sizes. For complex transforms, the recorded break-even is approximately:

- N=64: **431 transforms**;
- N=256: **131 transforms**;
- N=1024: **14 transforms**;
- N=4096: **1.4 transforms**;
- N≥16384: FFTW `ESTIMATE` is already cheaper to set up and faster to execute in this run.

Therefore a one-shot or short-run benchmark can have a different winner from a throughput benchmark. `FFTW_MEASURE` being the fastest execution kernel does not imply it minimizes end-to-end workload cost.

### Methodological conclusion

A defensible production-library comparison must specify at least:

- transform size and real/complex representation;
- planning policy;
- setup reuse count;
- normalization convention;
- in-place/out-of-place and allocation semantics;
- alignment/workspace rules;
- precision and thread count;
- exact library build and hardware.

See `pr4-fftw-baseline/ANALYSIS.md` and `docs/VENDOR_BENCHMARKS.md`.

## v3 — reusable planning and real-input specialization

The v3 corpus contains **2,790 raw observations**. It established that the same radix-2 decomposition becomes **1.295–1.593× faster** when reusable permutation/twiddle setup is removed from steady-state execution, with plan construction amortizing after roughly **1.75–3.95 transforms**. The specialized real path reaches roughly **1.74–1.75×** the planned complex throughput at the larger tested sizes.

The key conclusion was methodological: setup, complex execution, and real execution are separate workloads.

## v2 — algorithm families, structural complexity, and accuracy

The v2 research baseline established a broader algorithm taxonomy and independent accuracy studies. Selected results include:

- Rader about **1.29× faster than Bluestein** at N=509 and N=1009 in the recorded environment;
- iterative radix-2 about **1089× faster than direct DFT** at N=1024;
- a pedagogical split-radix implementation about **4.51× slower** than iterative radix-2 at N=256 despite a lower structural multiplication model;
- substantial algorithm-dependent forward-error differences across five deterministic signal families.

That milestone demonstrated that mathematical operation count, numerical behavior, and machine latency are related but not interchangeable.

## Evidence scope

All current formal baseline numbers come from virtualized/containerized development environments. They support implementation- and methodology-specific conclusions, not universal hardware rankings. The next external-validity step is repetition on named physical x86-64 and Arm machines with additional production backends and controlled threading/alignment.

The v1 baseline remains under `baseline-linux-amd-epyc-gcc14.csv` for historical continuity.
