# fft

A dependency-free **C++23 Fourier-transform laboratory**: multiple exact DFT/FFT algorithms implemented from first principles and compared on both **execution speed** and **numerical accuracy**.

The project is intentionally centered on one source file, [`fft.cpp`](fft.cpp), so the algorithms remain easy to inspect. It is a research/portfolio artifact rather than a wrapper around FFTW, MKL, Accelerate, cuFFT, or another tuned library.

## Algorithms

| Algorithm | Supported N | Time | Research role |
|---|---|---:|---|
| Direct DFT | any | O(N²) | mathematical baseline / oracle cross-check |
| Iterative radix-2 Cooley-Tukey | powers of two | O(N log N) | primary power-of-two implementation |
| Recursive radix-2 Cooley-Tukey | powers of two | O(N log N) | recursion/allocation comparison |
| Mixed-radix Cooley-Tukey | composite; prime leaves allowed | typically O(N log N) | factorization experiments |
| Bluestein / chirp-z | any, including primes | O(M log M) | efficient arbitrary-length transform |
| Auto dispatcher | any | varies | inspectable practical policy |

Forward and normalized inverse transforms are supported. The executable also contains deterministic signal generation, a long-double direct-DFT oracle, round-trip verification, an adaptive benchmark harness, and integrated self-tests.

## Build

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
./build/fft --self-test
```

Requires CMake 3.20+ and a C++23 compiler.

For host-specific benchmark builds:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DFFT_NATIVE=ON
cmake --build build -j
```

`FFT_NATIVE=ON` is useful for measuring one machine, but makes results less portable.

## Verify accuracy

```bash
./build/fft --verify --algorithm radix2-iterative --size 1024
./build/fft --verify --algorithm bluestein --size 1009
./build/fft --verify --algorithm mixed-radix --size 1000
```

Verification compares the double-precision result with an independently computed **long-double O(N²) DFT** and reports maximum absolute error, RMS error, relative L2 error, and forward+inverse round-trip error. The CLI caps oracle verification at N=4096 because it is deliberately quadratic.

## Benchmark

One transform:

```bash
./build/fft --benchmark --algorithm radix2-iterative --size 4096
./build/fft --benchmark --algorithm bluestein --size 4093 --samples 51 --target-ms 10
```

Full comparison:

```bash
./build/fft --benchmark-suite
./build/fft --benchmark-suite --csv > results.csv
```

The default harness performs untimed warmups, adaptively chooses repetitions so timing samples are long enough to measure, then collects **31 independent samples**. It reports min, p05, median, mean, p95, max, and sample standard deviation. A volatile checksum consumes transform output to prevent dead-code elimination.

Median is the primary latency statistic; p05/p95 and standard deviation make scheduler noise visible instead of hiding it behind one best run.

## Empirical baseline

[`results/baseline-linux-amd-epyc-gcc14.csv`](results/baseline-linux-amd-epyc-gcc14.csv) contains a 31-sample development baseline. [`results/SUMMARY.md`](results/SUMMARY.md) summarizes the main findings.

The baseline demonstrates several important effects immediately:

- iterative radix-2 decisively beats the direct DFT even at small powers of two in this implementation;
- recursive radix-2 is consistently slower because of allocations and recursion overhead;
- at prime N=1009, Bluestein is roughly two orders of magnitude faster than direct DFT and the mixed-radix implementation's prime leaf;
- the performance gap between O(N²) DFT and O(N log N) radix-2 rapidly becomes enormous as N grows.

These are **environment-specific measurements**, not universal rankings. Re-run on each target CPU/compiler before changing dispatch decisions.

## Auto policy

`auto` is intentionally simple:

- power of two → iterative radix-2;
- 2/3/5-smooth length → mixed-radix;
- otherwise → Bluestein.

This policy is exposed as a research baseline, not claimed to be globally optimal. A future milestone can learn crossover thresholds empirically per architecture/compiler.

## Reproducible benchmark reporting

When publishing results, record at least CPU model/architecture, OS, compiler/version, optimization flags, `FFT_NATIVE`, power mode, sample count, target sample duration, sizes, and algorithms. For cleaner experiments, use an idle machine, avoid thermal throttling, and repeat the complete suite across multiple sessions.

For deeper performance research, pair this harness with `perf`, Instruments, VTune, or equivalent tools and measure cycles, instructions, cache misses, and branch behavior. Platform counters are intentionally not embedded in this dependency-free core.

## Validation

```bash
ctest --test-dir build --output-on-failure
```

Sanitizers on GCC/Clang:

```bash
cmake -S . -B build-asan -DCMAKE_BUILD_TYPE=Debug -DFFT_SANITIZERS=ON
cmake --build build-asan -j
ctest --test-dir build-asan --output-on-failure
```

The self-test cross-checks arbitrary, prime, composite, and power-of-two sizes against direct DFT results, tests forward/inverse round trips, and verifies invalid radix-2 inputs are rejected.

## Scope

This is an algorithm laboratory, not a production replacement for FFTW, Intel oneMKL, Apple Accelerate/vDSP, cuFFT, rocFFT, or pocketfft. Those libraries are valuable future benchmark targets; fair comparisons must control planning time, allocation, threading, normalization, precision, and in-place/out-of-place behavior.

## License

MIT. See [LICENSE](LICENSE).
