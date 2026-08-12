# Benchmark summary

Development baseline: GCC 14.2.0, CMake Release, `FFT_NATIVE=OFF`, 31 timing samples per algorithm/size.

| N | Fastest implementation | Median | Runner-up | DFT speedup | p95 / median |
|---:|---|---:|---|---:|---:|
| 8 | `radix2-iterative` | 76.8 ns | radix2-recursive (284.1 ns) | 11.2× | 1.020 |
| 16 | `radix2-iterative` | 147.8 ns | radix2-recursive (615.3 ns) | 24.9× | 1.113 |
| 32 | `radix2-iterative` | 304.3 ns | radix2-recursive (1.31 µs) | 50.8× | 1.022 |
| 64 | `radix2-iterative` | 648.6 ns | radix2-recursive (2.74 µs) | 96.6× | 1.022 |
| 128 | `radix2-iterative` | 1.43 µs | radix2-recursive (5.81 µs) | 174.9× | 1.029 |
| 256 | `radix2-iterative` | 3.14 µs | radix2-recursive (12.23 µs) | 347.8× | 1.035 |
| 512 | `radix2-iterative` | 6.99 µs | radix2-recursive (25.51 µs) | 620.9× | 1.072 |
| 1009 | `bluestein` | 160.90 µs | direct DFT (16.55 ms) | 102.9× | 1.036 |
| 1024 | `radix2-iterative` | 15.48 µs | radix2-recursive (53.43 µs) | 1106.7× | 1.049 |

Across the recorded sizes, the geometric-mean speedup of the fastest measured implementation over direct DFT is about **121.8×**.

The results demonstrate the expected complexity crossover, the cost of recursive allocation-heavy radix-2, and the value of Bluestein for prime lengths. These figures are descriptive for this environment only; re-run the benchmark on each target machine before changing dispatch thresholds.
