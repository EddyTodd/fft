# Arbitrary-length plan analysis

Raw observations: **4,560**. Bootstrap seed: `20260812`; repetitions: **5,000**.

## Prime-length execution

| N | Bluestein M | Rader M | Legacy Bluestein | Planned Bluestein | Legacy Rader | Planned Rader | Planned winner | FFTW ESTIMATE | FFTW MEASURE |
|---:|---:|---:|---:|---:|---:|---:|---|---:|---:|
| 17 | 64 | 16 (cyclic) | 2.66 us | 744.9 ns | 1.31 us | 171.6 ns | **Rader** | 78.0 ns | 79.3 ns |
| 31 | 64 | 64 | 3.43 us | 767.1 ns | 2.55 us | 753.2 ns | **Rader** | 230.5 ns | 231.3 ns |
| 61 | 128 | 128 | 7.10 us | 1.65 us | 5.26 us | 1.66 us | **Bluestein** | 731.2 ns | 685.2 ns |
| 127 | 256 | 256 | 14.96 us | 3.63 us | 11.04 us | 3.68 us | **Bluestein** | 1.38 us | 1.13 us |
| 257 | 1024 | 256 (cyclic) | 52.43 us | 16.96 us | 23.72 us | 3.92 us | **Rader** | 3.18 us | 3.07 us |
| 509 | 1024 | 1024 | 66.75 us | 17.39 us | 49.48 us | 17.67 us | **Bluestein** | 6.09 us | 5.08 us |
| 1009 | 2048 | 2048 | 142.14 us | 38.02 us | 110.21 us | 38.62 us | **Bluestein** | 19.62 us | 19.68 us |
| 4093 | 8192 | 8192 | 680.10 us | 216.11 us | 512.46 us | 218.05 us | **unresolved** | 110.49 us | 59.17 us |

### Pairwise interpretation

- **N=17:** planned Bluestein / planned Rader median ratio 4.340x, 95% CI [4.299, 4.377] -> Rader. Planning improves Bluestein 3.57x and Rader 7.61x over the legacy setup-inclusive APIs. Best planned fftlab remains 2.16x slower than FFTW MEASURE.
- **N=31:** planned Bluestein / planned Rader median ratio 1.018x, 95% CI [1.008, 1.027] -> Rader. Planning improves Bluestein 4.47x and Rader 3.39x over the legacy setup-inclusive APIs. Best planned fftlab remains 3.26x slower than FFTW MEASURE.
- **N=61:** planned Bluestein / planned Rader median ratio 0.992x, 95% CI [0.984, 0.998] -> Bluestein. Planning improves Bluestein 4.31x and Rader 3.17x over the legacy setup-inclusive APIs. Best planned fftlab remains 2.40x slower than FFTW MEASURE.
- **N=127:** planned Bluestein / planned Rader median ratio 0.986x, 95% CI [0.980, 0.995] -> Bluestein. Planning improves Bluestein 4.12x and Rader 3.00x over the legacy setup-inclusive APIs. Best planned fftlab remains 3.20x slower than FFTW MEASURE.
- **N=257:** planned Bluestein / planned Rader median ratio 4.327x, 95% CI [4.278, 4.369] -> Rader. Planning improves Bluestein 3.09x and Rader 6.05x over the legacy setup-inclusive APIs. Best planned fftlab remains 1.28x slower than FFTW MEASURE.
- **N=509:** planned Bluestein / planned Rader median ratio 0.984x, 95% CI [0.967, 0.999] -> Bluestein. Planning improves Bluestein 3.84x and Rader 2.80x over the legacy setup-inclusive APIs. Best planned fftlab remains 3.43x slower than FFTW MEASURE.
- **N=1009:** planned Bluestein / planned Rader median ratio 0.984x, 95% CI [0.965, 0.997] -> Bluestein. Planning improves Bluestein 3.74x and Rader 2.85x over the legacy setup-inclusive APIs. Best planned fftlab remains 1.93x slower than FFTW MEASURE.
- **N=4093:** planned Bluestein / planned Rader median ratio 0.991x, 95% CI [0.976, 1.013] -> unresolved. Planning improves Bluestein 3.15x and Rader 2.35x over the legacy setup-inclusive APIs. Best planned fftlab remains 3.65x slower than FFTW MEASURE.

## Planning cost and amortization

| N | Bluestein setup | Bluestein break-even vs legacy | Rader setup | Rader break-even vs legacy | FFTW ESTIMATE setup | FFTW MEASURE setup |
|---:|---:|---:|---:|---:|---:|---:|
| 17 | 3.62 us | 1.89 | 3.75 us | 3.30 | 91.60 us | 524.11 us |
| 31 | 5.63 us | 2.11 | 4.14 us | 2.30 | 205.43 us | 67.563 ms |
| 61 | 6.67 us | 1.22 | 6.00 us | 1.67 | 855.05 us | 187.216 ms |
| 127 | 16.46 us | 1.45 | 11.45 us | 1.55 | 348.25 us | 153.613 ms |
| 257 | 32.86 us | 0.93 | 15.72 us | 0.79 | 730.79 us | 579.347 ms |
| 509 | 45.30 us | 0.92 | 34.23 us | 1.08 | 268.15 us | 283.139 ms |
| 1009 | 80.12 us | 0.77 | 84.81 us | 1.18 | 250.24 us | 186.670 ms |
| 4093 | 361.03 us | 0.78 | 288.56 us | 0.98 | 793.15 us | 804.286 ms |

## Interpretation

- Bluestein and Rader are both reduced to power-of-two convolution kernels so their planned comparison isolates reduction structure and persistent precomputation more cleanly than the legacy allocation-heavy APIs.
- When `N-1` is a power of two, planned Rader performs the cyclic convolution directly at length `N-1`; otherwise it zero-pads a linear convolution and folds it modulo `N-1`. That structural distinction is reported explicitly in the table.
- A lower convolution length is expected to help, but equal convolution lengths do not imply equal machine performance because Rader adds permutations/folding while Bluestein adds chirp multiplications.
- Setup and execution are separate response variables. The legacy functions intentionally remain setup-inclusive historical API baselines; planned-vs-planned and planned-vs-FFTW comparisons use persistent plans.
- FFTW is an adaptive production library and remains a reference point, not a claim that the local reductions have equivalent engineering maturity.
- Results are scoped to the recorded compiler/runtime/virtualized host and should not be converted into a universal prime-size dispatcher without physical-hardware replication.
