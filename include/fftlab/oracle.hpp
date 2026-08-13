#pragma once

#include "fftlab/fft.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <numbers>
#include <span>
#include <stdexcept>
#include <vector>

namespace fftlab {

using OracleComplex = std::complex<long double>;
using OracleVector = std::vector<OracleComplex>;

struct ErrorNorms {
    long double l1{};
    long double l2{};
    long double linf{};
};

template <FftScalar T>
[[nodiscard]] inline OracleVector oracle_dft(std::span<const ComplexT<T>> input,
                                             Direction direction = Direction::Forward) {
    const auto n = input.size();
    OracleVector output(n);
    if (n == 0) return output;
    const long double sign = direction == Direction::Forward ? -1.0L : 1.0L;
    const long double tau = 2.0L * std::numbers::pi_v<long double>;
    for (std::size_t k = 0; k < n; ++k) {
        OracleComplex sum{};
        for (std::size_t t = 0; t < n; ++t) {
            const auto angle = sign * tau * static_cast<long double>(k) *
                               static_cast<long double>(t) / static_cast<long double>(n);
            const OracleComplex x{static_cast<long double>(input[t].real()),
                                  static_cast<long double>(input[t].imag())};
            sum += x * OracleComplex{std::cos(angle), std::sin(angle)};
        }
        output[k] = sum;
    }
    if (direction == Direction::Inverse) {
        const auto scale = 1.0L / static_cast<long double>(n);
        for (auto& z : output) z *= scale;
    }
    return output;
}

template <FftScalar T>
[[nodiscard]] inline OracleVector widen(std::span<const ComplexT<T>> input) {
    OracleVector output;
    output.reserve(input.size());
    for (const auto z : input)
        output.emplace_back(static_cast<long double>(z.real()), static_cast<long double>(z.imag()));
    return output;
}

[[nodiscard]] inline ErrorNorms error_norms(std::span<const OracleComplex> actual,
                                             std::span<const OracleComplex> reference) {
    if (actual.size() != reference.size()) throw std::invalid_argument("error vector size mismatch");
    long double e1 = 0.0L, e2 = 0.0L, ei = 0.0L;
    long double r1 = 0.0L, r2 = 0.0L, ri = 0.0L;
    for (std::size_t i = 0; i < actual.size(); ++i) {
        const auto e = std::abs(actual[i] - reference[i]);
        const auto r = std::abs(reference[i]);
        e1 += e; e2 += e * e; ei = std::max(ei, e);
        r1 += r; r2 += r * r; ri = std::max(ri, r);
    }
    return {
        r1 != 0.0L ? e1 / r1 : e1,
        r2 != 0.0L ? std::sqrt(e2 / r2) : std::sqrt(e2),
        ri != 0.0L ? ei / ri : ei,
    };
}

template <FftScalar T>
[[nodiscard]] inline ErrorNorms error_norms(std::span<const ComplexT<T>> actual,
                                             std::span<const OracleComplex> reference) {
    const auto widened = widen<T>(actual);
    return error_norms(widened, reference);
}

} // namespace fftlab
