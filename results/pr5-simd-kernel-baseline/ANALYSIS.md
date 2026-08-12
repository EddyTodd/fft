# SIMD kernel analysis

Raw observations: **4,032**. Bootstrap seed: `20260812`; repetitions: **5,000**.

## Steady-state execution

| N | v3 plan | Scalar codelet | AVX2/FMA | AVX-512/FMA | Auto | FFTW ESTIMATE | FFTW MEASURE | Best SIMD / v3 | Gap to MEASURE closed |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 64 | 313.6 ns | 294.4 ns | 179.4 ns | 215.0 ns | 179.4 ns | 136.6 ns | 82.1 ns | **1.747x** (avx2) | **57.9%** |
| 256 | 1.60 us | 1.53 us | 758.1 ns | 936.4 ns | 759.2 ns | 464.0 ns | 367.5 ns | **2.104x** (avx2) | **68.2%** |
| 1024 | 7.62 us | 7.30 us | 3.75 us | 4.41 us | 3.75 us | 2.15 us | 1.70 us | **2.035x** (avx2) | **65.4%** |
| 4096 | 41.13 us | 38.22 us | 24.46 us | 21.07 us | 25.52 us | 15.08 us | 10.69 us | **1.952x** (avx512) | **65.9%** |
| 16384 | 228.53 us | 202.27 us | 142.49 us | 121.85 us | 142.80 us | 110.10 us | 51.76 us | **1.876x** (avx512) | **60.4%** |
| 65536 | 1.227 ms | 993.61 us | 733.51 us | 642.54 us | 719.95 us | 355.34 us | 292.70 us | **1.909x** (avx512) | **62.5%** |

### AVX2 vs AVX-512 and auto-selection

- **N=64:** AVX2/AVX-512 median ratio 0.835x, 95% CI [0.829, 0.843] -> AVX2 faster. Best SIMD speedup over v3 plan 1.747x, 95% CI [1.729, 1.759]. Best SIMD remains 2.19x slower than FFTW MEASURE. Auto selections: avx2: 3/3.
- **N=256:** AVX2/AVX-512 median ratio 0.810x, 95% CI [0.795, 0.854] -> AVX2 faster. Best SIMD speedup over v3 plan 2.104x, 95% CI [1.991, 2.142]. Best SIMD remains 2.06x slower than FFTW MEASURE. Auto selections: avx2: 3/3.
- **N=1024:** AVX2/AVX-512 median ratio 0.849x, 95% CI [0.842, 0.871] -> AVX2 faster. Best SIMD speedup over v3 plan 2.035x, 95% CI [1.992, 2.056]. Best SIMD remains 2.20x slower than FFTW MEASURE. Auto selections: avx2: 2/3, avx512: 1/3.
- **N=4096:** AVX2/AVX-512 median ratio 1.161x, 95% CI [1.157, 1.176] -> AVX-512 faster. Best SIMD speedup over v3 plan 1.952x, 95% CI [1.934, 1.968]. Best SIMD remains 1.97x slower than FFTW MEASURE. Auto selections: avx2: 2/3, avx512: 1/3.
- **N=16384:** AVX2/AVX-512 median ratio 1.169x, 95% CI [1.157, 1.178] -> AVX-512 faster. Best SIMD speedup over v3 plan 1.876x, 95% CI [1.852, 1.892]. Best SIMD remains 2.35x slower than FFTW MEASURE. Auto selections: avx2: 2/3, avx512: 1/3.
- **N=65536:** AVX2/AVX-512 median ratio 1.142x, 95% CI [1.129, 1.158] -> AVX-512 faster. Best SIMD speedup over v3 plan 1.909x, 95% CI [1.887, 1.938]. Best SIMD remains 2.20x slower than FFTW MEASURE. Auto selections: avx512: 3/3.

## Planning cost

| N | v3 setup | Best explicit SIMD setup | SIMD break-even vs v3 | Auto-tuned setup | Auto premium vs best explicit | FFTW ESTIMATE setup | FFTW MEASURE setup |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 64 | 3.83 us | 6.55 us (avx2) | 20.31 transforms | 26.483 ms | 26.476 ms | 132.82 us | 50.145 ms |
| 256 | 3.68 us | 11.68 us (avx2) | 9.56 transforms | 34.349 ms | 34.337 ms | 198.41 us | 137.230 ms |
| 1024 | 14.75 us | 24.05 us (avx2) | 2.40 transforms | 42.312 ms | 42.288 ms | 141.92 us | 249.103 ms |
| 4096 | 52.25 us | 79.36 us (avx512) | 1.35 transforms | 59.101 ms | 59.021 ms | 133.30 us | 483.840 ms |
| 16384 | 190.75 us | 289.80 us (avx512) | 0.93 transforms | 80.387 ms | 80.097 ms | 141.02 us | 1157.832 ms |
| 65536 | 1.046 ms | 1.436 ms (avx512) | 0.67 transforms | 98.930 ms | 97.493 ms | 770.57 us | 1869.938 ms |

## Interpretation

- The scalar codelet, AVX2, and AVX-512 paths share the same swap-list permutation and stage-contiguous twiddle layout. Their differences isolate instruction-width/code-generation effects more cleanly than comparing unrelated FFT algorithms.
- `KernelRadix2Plan` stores approximately N-1 stage-local complex twiddles rather than N/2 globally indexed twiddles. The additional persistent memory is part of the setup tradeoff and is not hidden from the analysis.
- The auto policy performs five rotated timing rounds over every supported candidate during plan construction. Auto tuning is therefore a portability/planning feature, not free execution speed; its full construction cost is reported separately.
- AVX2 and AVX-512 are explicit research modes even when auto chooses one of them. Wider vectors are only called faster when the bootstrap speedup interval excludes parity.
- FFTW remains a separately engineered adaptive library. Closing part of the latency gap with one vectorized radix-2 codelet does not imply architectural equivalence or a universal ranking.
- All findings are scoped to the recorded virtualized environment and compiler/runtime configuration.
