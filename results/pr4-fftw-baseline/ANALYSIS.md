# FFTW comparison analysis

Raw observations: **3,456**. Bootstrap seed: `20260812`.

## Complex execution

| N | fftlab planned | FFTW ESTIMATE | ESTIMATE speedup (95% CI) | FFTW MEASURE | MEASURE speedup (95% CI) | MEASURE / ESTIMATE |
|---:|---:|---:|---:|---:|---:|---:|
| 64 | 313.7 ns | 137.3 ns | **2.29×** [2.26, 2.31] | 81.6 ns | **3.84×** [3.82, 3.87] | **1.682×** [1.667, 1.695] |
| 256 | 1.58 µs | 468.3 ns | **3.38×** [3.23, 3.41] | 370.8 ns | **4.26×** [4.23, 4.32] | **1.263×** [1.250, 1.321] |
| 1024 | 7.77 µs | 2.13 µs | **3.66×** [3.63, 3.67] | 1.66 µs | **4.68×** [4.65, 4.70] | **1.280×** [1.274, 1.284] |
| 4096 | 40.55 µs | 15.01 µs | **2.70×** [2.68, 2.72] | 10.77 µs | **3.77×** [3.75, 3.78] | **1.394×** [1.385, 1.404] |
| 16384 | 230.84 µs | 108.18 µs | **2.13×** [2.12, 2.15] | 51.97 µs | **4.44×** [4.36, 4.50] | **2.082×** [2.049, 2.110] |
| 65536 | 1.20 ms | 354.77 µs | **3.39×** [3.34, 3.43] | 287.35 µs | **4.18×** [4.03, 4.29] | **1.235×** [1.191, 1.268] |

## Real execution

| N | fftlab planned | FFTW ESTIMATE | ESTIMATE speedup (95% CI) | FFTW MEASURE | MEASURE speedup (95% CI) | MEASURE / ESTIMATE |
|---:|---:|---:|---:|---:|---:|---:|
| 64 | 305.5 ns | 85.6 ns | **3.57×** [3.52, 3.61] | 66.8 ns | **4.57×** [4.52, 4.63] | **1.281×** [1.269, 1.295] |
| 256 | 1.35 µs | 323.6 ns | **4.17×** [4.13, 4.19] | 204.5 ns | **6.59×** [6.38, 6.63] | **1.582×** [1.531, 1.595] |
| 1024 | 6.03 µs | 1.03 µs | **5.85×** [5.83, 5.89] | 918.6 ns | **6.57×** [6.50, 6.64] | **1.122×** [1.110, 1.132] |
| 4096 | 27.35 µs | 5.13 µs | **5.34×** [5.30, 5.38] | 4.72 µs | **5.80×** [5.73, 5.84] | **1.086×** [1.074, 1.093] |
| 16384 | 138.47 µs | 41.48 µs | **3.34×** [3.32, 3.35] | 25.11 µs | **5.51×** [5.46, 5.56] | **1.652×** [1.635, 1.667] |
| 65536 | 711.37 µs | 211.74 µs | **3.36×** [3.26, 3.43] | 171.15 µs | **4.16×** [3.95, 4.28] | **1.237×** [1.176, 1.286] |

## Planning economics

| Kind | N | fftlab setup | FFTW ESTIMATE setup | ESTIMATE break-even vs fftlab | FFTW MEASURE cold setup | MEASURE break-even vs ESTIMATE |
|---|---:|---:|---:|---:|---:|---:|
| complex | 64 | 15.64 µs | 91.62 µs | 430.6 transforms | 49.36 ms | 885,401 transforms |
| complex | 256 | 7.73 µs | 152.94 µs | 130.5 transforms | 134.54 ms | 1,379,287 transforms |
| complex | 1024 | 14.15 µs | 91.29 µs | 13.7 transforms | 251.93 ms | 541,190 transforms |
| complex | 4096 | 59.38 µs | 94.17 µs | 1.4 transforms | 486.53 ms | 114,699 transforms |
| complex | 16384 | 205.06 µs | 133.40 µs | immediate | 1105.68 ms | 19,668 transforms |
| complex | 65536 | 949.91 µs | 737.70 µs | immediate | 1868.32 ms | 27,701 transforms |
| real | 64 | 2.34 µs | 193.49 µs | 869.0 transforms | 21.18 ms | 1,117,174 transforms |
| real | 256 | 18.51 µs | 556.43 µs | 524.8 transforms | 105.78 ms | 883,593 transforms |
| real | 1024 | 23.19 µs | 1.11 ms | 216.5 transforms | 239.43 ms | 2,126,344 transforms |
| real | 4096 | 46.41 µs | 1.53 ms | 66.7 transforms | 406.69 ms | 993,168 transforms |
| real | 16384 | 188.98 µs | 1.89 ms | 17.5 transforms | 782.12 ms | 47,651 transforms |
| real | 65536 | 790.61 µs | 2.36 ms | 3.1 transforms | 2048.21 ms | 50,403 transforms |

## Amortized winner

Winner minimizes `setup + K × execution` under the normalized forward/inverse contract.

| Kind | N | K=1 | K=10 | K=100 | K=10,000 | K=1,000,000 |
|---|---:|---|---|---|---|---|
| complex | 64 | fftlab | fftlab | fftlab | FFTW_ESTIMATE | FFTW_MEASURE |
| complex | 256 | fftlab | fftlab | fftlab | FFTW_ESTIMATE | FFTW_ESTIMATE |
| complex | 1024 | fftlab | fftlab | FFTW_ESTIMATE | FFTW_ESTIMATE | FFTW_MEASURE |
| complex | 4096 | fftlab | FFTW_ESTIMATE | FFTW_ESTIMATE | FFTW_ESTIMATE | FFTW_MEASURE |
| complex | 16384 | FFTW_ESTIMATE | FFTW_ESTIMATE | FFTW_ESTIMATE | FFTW_ESTIMATE | FFTW_MEASURE |
| complex | 65536 | FFTW_ESTIMATE | FFTW_ESTIMATE | FFTW_ESTIMATE | FFTW_ESTIMATE | FFTW_MEASURE |
| real | 64 | fftlab | fftlab | fftlab | FFTW_ESTIMATE | FFTW_ESTIMATE |
| real | 256 | fftlab | fftlab | fftlab | FFTW_ESTIMATE | FFTW_MEASURE |
| real | 1024 | fftlab | fftlab | fftlab | FFTW_ESTIMATE | FFTW_ESTIMATE |
| real | 4096 | fftlab | fftlab | FFTW_ESTIMATE | FFTW_ESTIMATE | FFTW_MEASURE |
| real | 16384 | fftlab | fftlab | FFTW_ESTIMATE | FFTW_ESTIMATE | FFTW_MEASURE |
| real | 65536 | fftlab | FFTW_ESTIMATE | FFTW_ESTIMATE | FFTW_ESTIMATE | FFTW_MEASURE |

## Effect-size checks

Common-language values below are the probability that a random FFTW MEASURE execution observation is faster than a random fftlab planned observation.

| Kind | N | P(FFTW MEASURE faster) |
|---|---:|---:|
| complex | 64 | 1.0000 |
| complex | 256 | 0.9892 |
| complex | 1024 | 1.0000 |
| complex | 4096 | 1.0000 |
| complex | 16384 | 1.0000 |
| complex | 65536 | 0.9896 |
| real | 64 | 1.0000 |
| real | 256 | 1.0000 |
| real | 1024 | 1.0000 |
| real | 4096 | 1.0000 |
| real | 16384 | 1.0000 |
| real | 65536 | 1.0000 |

## Interpretation

- FFTW execution is normalized to the same mathematical contract as fftlab: forward plus normalized inverse. FFTW’s required `1/N` inverse scaling is included in its timed path.
- FFTW planning arrays are allocated before setup timing. `FFTW_MEASURE` setup is cold: wisdom is forgotten before each measured forward+inverse plan pair.
- `FFTW_ESTIMATE` and `FFTW_MEASURE` execution use persistent plans and preallocated buffers; planning is excluded from execution timing.
- A MEASURE execution win is not automatically an end-to-end win. The break-even table quantifies how many repeated transforms are required to repay additional planning time.
- Results describe this binary/library/virtualized machine. They do not imply FFTW will lead by the same factor on another architecture, library build, compiler, transform shape, thread count, or planning policy.
