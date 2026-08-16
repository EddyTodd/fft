#include "fftlab/oracle.hpp"
#include "fftlab/planner.hpp"

#include <cmath>
#include <iostream>
#include <random>
#include <stdexcept>
#include <vector>

namespace {
template <fftlab::FftScalar T>
bool close(fftlab::ComplexT<T> a, fftlab::ComplexT<T> b) {
    const T eps = std::same_as<T, float> ? T{2e-3} : T{2e-10};
    return std::abs(a - b) <= eps * (T{1} + std::abs(b));
}

void require(bool ok, const char* message) { if (!ok) throw std::runtime_error(message); }

template <fftlab::FftScalar T>
void exercise(std::size_t n, fftlab::PlanPreference preference) {
    fftlab::Plan<T> plan(n, {preference, true});
    std::mt19937_64 rng{0x504c414eULL + n};
    std::uniform_real_distribution<double> dist{-1.0, 1.0};
    fftlab::VectorT<T> input(n);
    for (auto& z : input) z = {static_cast<T>(dist(rng)), static_cast<T>(dist(rng))};

    const auto oracle = fftlab::oracle_dft<T>(input);
    fftlab::VectorT<T> output(n), recovered(n), scratch(plan.scratch_size());
    plan.forward(input, output, scratch);
    for (std::size_t i = 0; i < n; ++i) {
        const fftlab::ComplexT<T> expected{static_cast<T>(oracle[i].real()), static_cast<T>(oracle[i].imag())};
        require(close(output[i], expected), "out-of-place plan disagrees with oracle");
    }
    plan.inverse(output, recovered, scratch);
    for (std::size_t i = 0; i < n; ++i) require(close(recovered[i], input[i]), "out-of-place roundtrip failed");

    auto inplace = input;
    fftlab::VectorT<T> inplace_scratch(plan.inplace_scratch_size());
    plan.forward_inplace(inplace, inplace_scratch);
    for (std::size_t i = 0; i < n; ++i) require(close(inplace[i], output[i]), "in-place forward diverged");
    plan.inverse_inplace(inplace, inplace_scratch);
    for (std::size_t i = 0; i < n; ++i) require(close(inplace[i], input[i]), "in-place roundtrip failed");

    if (plan.scratch_size() > 0) {
        bool rejected = false;
        try {
            fftlab::VectorT<T> short_scratch(plan.scratch_size() - 1);
            plan.forward(input, output, short_scratch);
        } catch (const std::invalid_argument&) { rejected = true; }
        require(rejected, "out-of-place plan accepted insufficient scratch");
    }
    if (plan.inplace_scratch_size() > 0) {
        bool rejected = false;
        try {
            fftlab::VectorT<T> short_scratch(plan.inplace_scratch_size() - 1);
            auto data = input;
            plan.forward_inplace(data, short_scratch);
        } catch (const std::invalid_argument&) { rejected = true; }
        require(rejected, "in-place plan accepted insufficient scratch");
    }

    bool rejected_size = false;
    try {
        fftlab::VectorT<T> wrong(n + 1);
        plan.forward_inplace(wrong, inplace_scratch);
    } catch (const std::invalid_argument&) { rejected_size = true; }
    require(rejected_size, "plan accepted wrong data size");
}

template <fftlab::FftScalar T>
void run_type() {
    exercise<T>(8, fftlab::PlanPreference::Radix2);
    exercise<T>(12, fftlab::PlanPreference::MixedRadix);
    exercise<T>(15, fftlab::PlanPreference::GoodThomas);
    exercise<T>(17, fftlab::PlanPreference::Rader);
    exercise<T>(31, fftlab::PlanPreference::Bluestein);

    fftlab::Plan<T> identity0(0);
    fftlab::VectorT<T> empty;
    identity0.forward(empty, empty);
    identity0.forward_inplace(empty);

    fftlab::Plan<T> identity1(1);
    fftlab::VectorT<T> one{{T{3}, T{-2}}}, out(1);
    identity1.forward(one, out);
    require(out == one, "N=1 identity failed");
}
}

int main() {
    try {
        run_type<float>();
        run_type<double>();
        std::cout << "PASS: public Plan<T> contract matrix\n";
    } catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
