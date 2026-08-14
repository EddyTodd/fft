# Contributing

`fftlab` owns permanent Fourier-transform algorithms, reusable plans, planner semantics, numerical correctness/oracle utilities, explicit SIMD kernels, and theory.

Generic timing, campaign orchestration, vendor performance comparisons, statistical analysis, provenance, and benchmark result corpora belong in `EddyTodd/bench`.

## Algorithm and planner changes

Document the supported structural domain, transform convention, allocation/scratch behavior, and planner-selection implications. New mechanisms require deterministic forward/inverse correctness checks against the independent oracle on representative valid lengths.

## SIMD changes

Keep ISA selection explicit and runtime-safe. Unsupported explicit ISA requests must fail rather than silently execute a differently named path. Validate equivalence with the scalar implementation on capable hardware.

## Validation

Run normal, release, sanitizer, header-self-containment, package-consumer, and release-metadata checks. Package changes must remain relocatable.

Do not add benchmark runners, vendor timing frontends, statistical analyzers, or result-management infrastructure to this repository.
