#include "fftlab/arbitrary_plan.hpp"

#include "internal.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <type_traits>

namespace fftlab {
namespace {

std::size_t bluestein_workspace(std::size_t n) {
    if (n == 0) throw std::invalid_argument("BluesteinPlan requires N >= 1");
    if (n == 1) return 1;
    if (n > (std::numeric_limits<std::size_t>::max() / 2) + 1)
        throw std::length_error("BluesteinPlan workspace overflow");
    return next_pow2(2 * n - 1);
}

std::size_t rader_workspace(std::size_t n) {
    const auto l = n - 1;
    if (pow2(l)) return l;
    if (l > (std::numeric_limits<std::size_t>::max() / 2) + 1)
        throw std::length_error("RaderPlan workspace overflow");
    return next_pow2(2 * l - 1);
}

Complex forward_chirp(std::size_t k, std::size_t n) {
    const long double kd = static_cast<long double>(k);
    const long double period = 2.0L * static_cast<long double>(n);
    const long double phase = std::fmod(kd * kd, period) / static_cast<long double>(n);
    return root(-pi * static_cast<double>(phase));
}

} // namespace

BluesteinPlan::BluesteinPlan(std::size_t n)
    : n_(n), m_(bluestein_workspace(n)), convolution_plan_(m_), chirp_(n), kernel_spectrum_(m_) {
    for (std::size_t k = 0; k < n_; ++k) {
        chirp_[k] = forward_chirp(k, n_);
        const auto b = std::conj(chirp_[k]);
        kernel_spectrum_[k] = b;
        if (k != 0) kernel_spectrum_[m_ - k] = b;
    }
    convolution_plan_.forward_inplace(kernel_spectrum_);
}

void BluesteinPlan::forward(const Vector& input, Vector& output, Vector& scratch) const {
    execute(input, output, scratch, false);
}

void BluesteinPlan::inverse(const Vector& input, Vector& output, Vector& scratch) const {
    execute(input, output, scratch, true);
}

void BluesteinPlan::execute(const Vector& input, Vector& output, Vector& scratch, bool inverse) const {
    if (input.size() != n_ || output.size() != n_ || scratch.size() != m_)
        throw std::invalid_argument("BluesteinPlan buffer size mismatch");

    std::fill(scratch.begin(), scratch.end(), Complex{});
    for (std::size_t k = 0; k < n_; ++k) {
        const auto x = inverse ? std::conj(input[k]) : input[k];
        scratch[k] = x * chirp_[k];
    }
    convolution_plan_.forward_inplace(scratch);
    for (std::size_t k = 0; k < m_; ++k) scratch[k] *= kernel_spectrum_[k];
    convolution_plan_.inverse_inplace(scratch);

    if (!inverse) {
        for (std::size_t k = 0; k < n_; ++k) output[k] = scratch[k] * chirp_[k];
    } else {
        const double scale = 1.0 / static_cast<double>(n_);
        for (std::size_t k = 0; k < n_; ++k) output[k] = std::conj(scratch[k] * chirp_[k]) * scale;
    }
}

RaderPlan::RaderPlan(std::size_t prime_n)
    : n_(prime_n), l_(prime_n > 0 ? prime_n - 1 : 0), m_([&] {
          if (prime_n < 3 || !is_prime(prime_n))
              throw std::invalid_argument("RaderPlan requires prime N >= 3");
          return rader_workspace(prime_n);
      }()),
      direct_cyclic_(m_ == l_), convolution_plan_(m_), output_permutation_(l_),
      input_permutation_(l_), kernel_spectrum_(m_) {
    const auto g = primitive_root_prime(n_);
    output_permutation_[0] = 1;
    for (std::size_t q = 1; q < l_; ++q)
        output_permutation_[q] = mul_mod(output_permutation_[q - 1], g, n_);
    for (std::size_t q = 0; q < l_; ++q)
        input_permutation_[q] = output_permutation_[(l_ - q) % l_];

    for (std::size_t q = 0; q < l_; ++q)
        kernel_spectrum_[q] = root(-2.0 * pi * static_cast<double>(output_permutation_[q]) /
                                   static_cast<double>(n_));
    convolution_plan_.forward_inplace(kernel_spectrum_);
}

void RaderPlan::forward(const Vector& input, Vector& output, Vector& scratch) const {
    execute(input, output, scratch, false);
}

void RaderPlan::inverse(const Vector& input, Vector& output, Vector& scratch) const {
    execute(input, output, scratch, true);
}

void RaderPlan::execute(const Vector& input, Vector& output, Vector& scratch, bool inverse) const {
    if (input.size() != n_ || output.size() != n_ || scratch.size() != m_)
        throw std::invalid_argument("RaderPlan buffer size mismatch");

    std::fill(scratch.begin(), scratch.end(), Complex{});
    Complex dc{};
    for (std::size_t i = 0; i < n_; ++i) {
        const auto x = inverse ? std::conj(input[i]) : input[i];
        dc += x;
    }
    for (std::size_t q = 0; q < l_; ++q) {
        const auto x = inverse ? std::conj(input[input_permutation_[q]])
                               : input[input_permutation_[q]];
        scratch[q] = x;
    }

    convolution_plan_.forward_inplace(scratch);
    for (std::size_t i = 0; i < m_; ++i) scratch[i] *= kernel_spectrum_[i];
    convolution_plan_.inverse_inplace(scratch);

    const auto x0 = inverse ? std::conj(input[0]) : input[0];
    if (!inverse) {
        output[0] = dc;
        for (std::size_t q = 0; q < l_; ++q) {
            auto c = scratch[q];
            if (!direct_cyclic_ && q + l_ < 2 * l_ - 1) c += scratch[q + l_];
            output[output_permutation_[q]] = x0 + c;
        }
    } else {
        const double scale = 1.0 / static_cast<double>(n_);
        output[0] = std::conj(dc) * scale;
        for (std::size_t q = 0; q < l_; ++q) {
            auto c = scratch[q];
            if (!direct_cyclic_ && q + l_ < 2 * l_ - 1) c += scratch[q + l_];
            output[output_permutation_[q]] = std::conj(x0 + c) * scale;
        }
    }
}

ArbitraryPlan::Storage ArbitraryPlan::make_storage(std::size_t n, ArbitraryPlanPolicy policy,
                                                    ArbitraryPlanAlgorithm& algorithm) {
    if (n == 0) throw std::invalid_argument("ArbitraryPlan requires N >= 1");
    if (policy == ArbitraryPlanPolicy::Rader) {
        algorithm = ArbitraryPlanAlgorithm::Rader;
        return Storage{std::in_place_type<RaderPlan>, n};
    }
    if (policy == ArbitraryPlanPolicy::Bluestein) {
        algorithm = ArbitraryPlanAlgorithm::Bluestein;
        return Storage{std::in_place_type<BluesteinPlan>, n};
    }
    if (pow2(n)) {
        algorithm = ArbitraryPlanAlgorithm::Radix2;
        return Storage{std::in_place_type<Radix2Plan>, n};
    }
    if (n >= 3 && is_prime(n)) {
        algorithm = ArbitraryPlanAlgorithm::Rader;
        return Storage{std::in_place_type<RaderPlan>, n};
    }
    algorithm = ArbitraryPlanAlgorithm::Bluestein;
    return Storage{std::in_place_type<BluesteinPlan>, n};
}

ArbitraryPlan::ArbitraryPlan(std::size_t n, ArbitraryPlanPolicy policy)
    : n_(n), storage_(make_storage(n, policy, algorithm_)) {}

std::size_t ArbitraryPlan::scratch_size() const noexcept {
    return std::visit([](const auto& p) -> std::size_t {
        using T = std::decay_t<decltype(p)>;
        if constexpr (std::is_same_v<T, Radix2Plan>) return 0;
        else return p.scratch_size();
    }, storage_);
}

void ArbitraryPlan::forward(const Vector& input, Vector& output, Vector& scratch) const {
    if (input.size() != n_ || output.size() != n_)
        throw std::invalid_argument("ArbitraryPlan buffer size mismatch");
    std::visit([&](const auto& p) {
        using T = std::decay_t<decltype(p)>;
        if constexpr (std::is_same_v<T, Radix2Plan>) {
            output = input;
            p.forward_inplace(output);
        } else {
            p.forward(input, output, scratch);
        }
    }, storage_);
}

void ArbitraryPlan::inverse(const Vector& input, Vector& output, Vector& scratch) const {
    if (input.size() != n_ || output.size() != n_)
        throw std::invalid_argument("ArbitraryPlan buffer size mismatch");
    std::visit([&](const auto& p) {
        using T = std::decay_t<decltype(p)>;
        if constexpr (std::is_same_v<T, Radix2Plan>) {
            output = input;
            p.inverse_inplace(output);
        } else {
            p.inverse(input, output, scratch);
        }
    }, storage_);
}

void arbitrary_plan_tests() {
    std::size_t checks = 0;
    auto req = [&](bool ok, const char* what) {
        ++checks;
        if (!ok) throw std::runtime_error(std::string("arbitrary plan self-test failure: ") + what);
    };
    std::mt19937_64 rng(0xA8B17A4ULL);
    std::uniform_real_distribution<double> dist(-1.0, 1.0);

    for (std::size_t n : {1u, 2u, 3u, 5u, 7u, 8u, 12u, 17u, 31u, 64u, 127u, 255u, 509u}) {
        Vector x(n);
        for (auto& z : x) z = {dist(rng), dist(rng)};
        const auto reference = dft(x);

        BluesteinPlan bp(n);
        Vector bout(n), bs(bp.scratch_size());
        bp.forward(x, bout, bs);
        for (std::size_t i = 0; i < n; ++i)
            req(std::abs(bout[i] - reference[i]) < 3e-10 * (1.0 + std::abs(reference[i])), "Bluestein forward");
        Vector back(n);
        bp.inverse(bout, back, bs);
        for (std::size_t i = 0; i < n; ++i)
            req(std::abs(back[i] - x[i]) < 3e-10 * (1.0 + std::abs(x[i])), "Bluestein roundtrip");

        ArbitraryPlan ap(n);
        Vector aout(n), as(ap.scratch_size());
        ap.forward(x, aout, as);
        for (std::size_t i = 0; i < n; ++i)
            req(std::abs(aout[i] - reference[i]) < 3e-10 * (1.0 + std::abs(reference[i])), "ArbitraryPlan forward");
        ap.inverse(aout, back, as);
        for (std::size_t i = 0; i < n; ++i)
            req(std::abs(back[i] - x[i]) < 3e-10 * (1.0 + std::abs(x[i])), "ArbitraryPlan roundtrip");

        if (n >= 3 && is_prime(n)) {
            RaderPlan rp(n);
            Vector rout(n), rs(rp.scratch_size());
            rp.forward(x, rout, rs);
            for (std::size_t i = 0; i < n; ++i)
                req(std::abs(rout[i] - reference[i]) < 3e-10 * (1.0 + std::abs(reference[i])), "Rader forward");
            rp.inverse(rout, back, rs);
            for (std::size_t i = 0; i < n; ++i)
                req(std::abs(back[i] - x[i]) < 3e-10 * (1.0 + std::abs(x[i])), "Rader roundtrip");
        }
    }

    bool threw = false;
    try { RaderPlan bad(15); } catch (const std::invalid_argument&) { threw = true; }
    req(threw, "Rader rejects composite N");
    threw = false;
    try { ArbitraryPlan bad(0); } catch (const std::invalid_argument&) { threw = true; }
    req(threw, "ArbitraryPlan rejects N=0");

    std::cout << "PASS: " << checks << " arbitrary-plan checks\n";
}

} // namespace fftlab
