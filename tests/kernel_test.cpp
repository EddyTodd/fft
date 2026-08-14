#include "fftlab/kernel.hpp"
#include "fftlab/plan.hpp"

#include <cmath>
#include <cstddef>
#include <iostream>
#include <random>
#include <stdexcept>
#include <vector>

using namespace fftlab;

namespace {

void require_close(Complex64 actual, Complex64 reference, const char* message,
                   std::size_t& checks) {
    ++checks;
    if (std::abs(actual - reference) > 3e-12 * (1.0 + std::abs(reference))) {
        throw std::runtime_error(message);
    }
}

}  // namespace

int main() {
    const auto capabilities = kernel_capabilities();
    std::size_t checks = 0;

    std::mt19937_64 rng{7};
    std::uniform_real_distribution<double> distribution{-1.0, 1.0};

    std::vector<KernelIsa> modes{KernelIsa::Scalar};
    if (capabilities.avx2) {
        modes.push_back(KernelIsa::Avx2);
    }
    if (capabilities.avx512) {
        modes.push_back(KernelIsa::Avx512);
    }

    for (std::size_t n : {1u, 2u, 4u, 8u, 16u, 64u, 256u, 1024u}) {
        Vector64 input(n);
        for (auto& value : input) {
            value = {distribution(rng), distribution(rng)};
        }

        Radix2Plan<double> reference_plan(n);
        auto reference = input;
        reference_plan.forward_inplace(reference);

        for (const auto isa : modes) {
            KernelRadix2Plan plan(n, isa);
            auto output = input;
            plan.forward_inplace(output);
            for (std::size_t i = 0; i < n; ++i) {
                require_close(output[i], reference[i], "kernel mismatch", checks);
            }

            plan.inverse_inplace(output);
            for (std::size_t i = 0; i < n; ++i) {
                require_close(output[i], input[i], "kernel roundtrip mismatch", checks);
            }
        }
    }

    std::cout << "PASS: " << checks << " kernel checks; avx2=" << capabilities.avx2
              << " avx512=" << capabilities.avx512 << '\n';
}
