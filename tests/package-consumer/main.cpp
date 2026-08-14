#include <fftlab/fftlab.hpp>

#include <cmath>
#include <numbers>
#include <span>
#include <vector>

int main() {
    constexpr std::size_t n = 16;
    fftlab::Vector64 input(n);
    for (std::size_t i = 0; i < n; ++i) {
        const double angle = 2.0 * std::numbers::pi_v<double> * static_cast<double>(i) /
                             static_cast<double>(n);
        input[i] = {std::cos(angle), std::sin(angle)};
    }

    fftlab::Plan<double> plan(n);
    fftlab::Vector64 output(n);
    fftlab::Vector64 scratch(plan.scratch_size());
    plan.forward(input, output, scratch);

    fftlab::Vector64 roundtrip(n);
    plan.inverse(output, roundtrip, scratch);
    for (std::size_t i = 0; i < n; ++i) {
        if (std::abs(roundtrip[i] - input[i]) > 1e-10 * (1.0 + std::abs(input[i]))) {
            return 1;
        }
    }

    return 0;
}
