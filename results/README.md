# Benchmark results

`baseline-linux-amd-epyc-gcc14.csv` is a development baseline captured while building v1. It makes the repository empirical from its first revision; it is not a claim of universal hardware rankings.

Environment:

- CPU model exposed to the container: AMD EPYC 9V74 80-Core Processor
- architecture: x86_64
- OS kernel: Linux 6.18.35
- compiler: GCC 14.2.0
- build: CMake `Release`, `FFT_NATIVE=OFF`
- samples: 31 per algorithm/size
- warmups: 5
- adaptive target duration: 1 ms per timing sample

The environment is virtualized/containerized, so absolute timing should not be treated as a physical-hardware reference. The CSV is primarily useful for validating the benchmark pipeline and preserving the observations that motivated the initial dispatcher.

For portfolio-quality hardware claims, regenerate the suite on a named physical machine using the methodology in the root README and commit its environment metadata beside the CSV.
