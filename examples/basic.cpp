#include <fftlab/fftlab.hpp>

#include <array>
#include <cmath>

int main() {
    const std::array<fftlab::Complex64, 4> input{{
        {1.0, 0.0},
        {2.0, 0.0},
        {3.0, 0.0},
        {4.0, 0.0},
    }};
    std::array<fftlab::Complex64, 4> spectrum{};
    std::array<fftlab::Complex64, 4> roundtrip{};

    fftlab::Plan<double> plan(input.size());
    fftlab::Vector64 scratch(plan.scratch_size());
    plan.forward(input, spectrum, scratch);
    plan.inverse(spectrum, roundtrip, scratch);

    for (std::size_t i = 0; i < input.size(); ++i) {
        if (std::abs(roundtrip[i] - input[i]) > 1e-12) {
            return 1;
        }
    }
    return 0;
}
