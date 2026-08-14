# v1 completeness and scope

## Definition of v1

fftlab 1.0 is complete when it provides a coherent representative catalog of important one-dimensional CPU DFT/FFT mechanisms, reusable structural planning for arbitrary lengths, binary32/binary64 support, deterministic numerical validation, and an installable C++23 library boundary independent of benchmarking infrastructure.

## Included

### Data and execution

- one-dimensional transforms;
- CPU execution;
- sequential algorithms/plans;
- complex binary32 and binary64;
- real power-of-two transforms in binary32 and binary64;
- arbitrary complex lengths;
- reusable plans and caller-owned scratch.

### Algorithm mechanisms

- direct DFT;
- iterative radix-2 Cooley-Tukey;
- recursive radix-2;
- Stockham autosort radix-2;
- radix-4;
- classical split-radix;
- scaled modified split-radix (Johnson-Frigo mechanism);
- mixed-radix Cooley-Tukey;
- Good-Thomas / Prime Factor Algorithm for coprime decompositions;
- Rader for prime lengths;
- Bluestein for arbitrary lengths;
- small planned DFT codelets for radices 2/3/4/5/7.

### Reusable planning

- radix-2 plans;
- real radix-2 plans;
- planned mixed-radix;
- Good-Thomas plans;
- Rader plans;
- Bluestein plans;
- arbitrary-length `Plan<T>` with explicit structural policy/capability queries;
- explicit binary64 scalar/AVX2/AVX-512 radix-2 kernel plans.

### Quality boundary

- C++23;
- strict warning-clean GCC/Clang builds;
- ASan/UBSan validation;
- deterministic long-double oracle suite for f32/f64;
- CMake install/export/package configuration;
- no required FFTW dependency for the installed core.

## Explicitly deferred after v1

The following are useful Fourier-transform domains, but they do not block the v1 1D CPU algorithm library:

- multidimensional transforms;
- batched transform APIs and layout planners;
- internal multithreading and distributed FFTs;
- GPU/CUDA/HIP/Metal FFTs;
- DCT/DST and other real-to-real transform families;
- nonuniform FFT (NUFFT);
- pruned/sparse FFT APIs;
- arbitrary-precision production transforms;
- SVE and dedicated Arm NEON codelets;
- dedicated binary32 x86 SIMD codelets;
- cache-oblivious/four-step/six-step large-transform schedulers;
- generated codelet systems;
- vendor backends as a core dependency;
- empirical auto-tuning/wisdom databases.

These are post-v1 feature domains, not unresolved correctness defects in the v1 contract.

## Why NEON is deferred

The portable f32/f64 algorithms and plans compile independently of x86 ISA support. A representative NEON kernel would be valuable, but it should be implemented and validated on actual Arm hardware rather than added untested from an x86-only environment. v1 therefore states this boundary explicitly. The scalar/generic planner remains the required Arm path.

## Benchmark/research separation

Development produced substantial benchmark executables, FFTW adapters, campaign scripts, statistical analyzers, and raw result sets. They are not part of the installed v1 package. They are historical/legacy research inputs whose long-term home is `EddyTodd/bench`.

See `LEGACY_RESEARCH.md` for the migration manifest.

## Core blocker status

At the v1 freeze, no known mechanism from the required representative catalog is intentionally left as documentation-only. Deferred items above expand domain, architecture specialization, or empirical infrastructure rather than filling a missing required v1 mechanism.
