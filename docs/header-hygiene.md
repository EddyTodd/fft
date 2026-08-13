# Public-header hygiene

`fftlab` treats every header in its 15-file `public_headers` CMake file set as independently includable C++23 API.

When `FFTLAB_BUILD_TESTS=ON`, CMake reads the target file set directly and generates one translation unit per header. Each generated source contains only that FFT header include. The `fftlab-header-self-containment` object target compiles all of those sources as part of the normal test-enabled build.

This catches accidental include-order coupling across algorithm, plan, planner, codelet, kernel, oracle, and type headers. Public headers must include every declaration they require rather than depending on a consumer to include `fftlab.hpp`, `types.hpp`, or another project header first.

The check complements `fftlab.package-consumer`: self-containment validates each header independently in the source tree, while the package consumer validates the complete installed/exported header surface from a separate CMake project.
