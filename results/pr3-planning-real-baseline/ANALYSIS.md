# Planned and real FFT analysis

Raw observations: **2,790** from **3** input file(s). Bootstrap seed: `20260812`.

| N | Legacy complex | Planned complex | Plan speedup (95% CI) | Planned real | Real speedup (95% CI) | Complex setup | Plan break-even | Real setup | Real break-even |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 64 | 661.3 ns | 415.2 ns | **1.593×** [1.578, 1.606] | 350.1 ns | **1.186×** [1.177, 1.196] | 971.0 ns | 3.95 transforms | 1.33 µs | 5.54 transforms |
| 256 | 3.27 µs | 2.12 µs | **1.542×** [1.522, 1.568] | 1.61 µs | **1.320×** [1.298, 1.339] | 3.15 µs | 2.74 transforms | 3.52 µs | 0.70 transforms |
| 1024 | 16.08 µs | 10.41 µs | **1.545×** [1.529, 1.582] | 7.35 µs | **1.415×** [1.368, 1.438] | 11.14 µs | 1.96 transforms | 11.80 µs | 0.22 transforms |
| 4096 | 77.65 µs | 51.65 µs | **1.504×** [1.482, 1.522] | 32.72 µs | **1.579×** [1.558, 1.602] | 45.41 µs | 1.75 transforms | 45.67 µs | 0.01 transforms |
| 16384 | 386.18 µs | 282.78 µs | **1.366×** [1.332, 1.406] | 161.27 µs | **1.753×** [1.701, 1.804] | 213.38 µs | 2.06 transforms | 209.59 µs | 0.00 transforms |
| 65536 | 1774.75 µs | 1370.09 µs | **1.295×** [1.273, 1.316] | 785.55 µs | **1.744×** [1.708, 1.773] | 1022.69 µs | 2.53 transforms | 937.85 µs | 0.00 transforms |

## Interpretation

- `legacy-complex` and `planned-complex` execute the same radix-2 decomposition; the planned path moves bit-reversal and trigonometric twiddle construction out of the timed kernel.
- `planned-real` exploits real-input Hermitian symmetry by packing even/odd samples into an N/2 complex FFT and reconstructing only N/2+1 unique frequency bins.
- Setup break-even is descriptive for this implementation and machine: plan construction cost divided by the median per-transform savings.
- The real-path break-even compares the extra real-plan setup cost against planned-complex setup; a value below one means the specialized real representation repays any extra setup within the first transform.
- Confidence intervals use independent nonparametric bootstrap resampling of the pooled raw observations. They characterize this recorded environment, not universal hardware rankings.
