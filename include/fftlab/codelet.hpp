#pragma once

#include "fftlab/plan.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <memory>
#include <span>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

namespace fftlab {

enum class PlanAlgorithm { Identity, Radix2, MixedRadix, GoodThomas, Rader, Bluestein };
enum class PlanPreference { Structural, Radix2, MixedRadix, GoodThomas, Rader, Bluestein };

struct PlanOptions {
    PlanPreference preference{PlanPreference::Structural};
    bool prefer_good_thomas{true};
};

struct PlanCapabilities {
    bool radix2{};
    bool mixed_radix{};
    bool good_thomas{};
    bool rader{};
    bool bluestein{true};
};

[[nodiscard]] inline bool mixed_radix_decomposable(std::size_t n) noexcept {
    if (n < 2) return false;
    for (const auto p : {2U, 3U, 5U, 7U}) while (n % p == 0) n /= p;
    return n == 1;
}

[[nodiscard]] inline PlanCapabilities plan_capabilities(std::size_t n) noexcept {
    const auto split = coprime_factor_split(n);
    return {
        n >= 1 && pow2(n),
        mixed_radix_decomposable(n) && !pow2(n),
        split.first > 1 && split.second > 1,
        n >= 3 && is_prime(n),
        n >= 1
    };
}

[[nodiscard]] inline std::string_view plan_name(PlanAlgorithm a) noexcept {
    switch (a) {
        case PlanAlgorithm::Identity: return "identity";
        case PlanAlgorithm::Radix2: return "radix2";
        case PlanAlgorithm::MixedRadix: return "mixed-radix";
        case PlanAlgorithm::GoodThomas: return "good-thomas";
        case PlanAlgorithm::Rader: return "rader";
        case PlanAlgorithm::Bluestein: return "bluestein";
    }
    return "unknown";
}

[[nodiscard]] inline constexpr bool small_codelet_supported(std::size_t radix) noexcept {
    return radix == 2 || radix == 3 || radix == 4 || radix == 5 || radix == 7;
}

[[nodiscard]] inline std::size_t choose_small_radix(std::size_t n) noexcept {
    // Radix-4 is intentionally preferred over two radix-2 levels; remaining
    // codelets cover the most useful small prime factors for CPU FFT planning.
    for (const auto r : {4U, 2U, 3U, 5U, 7U}) if (n % r == 0) return r;
    return 0;
}

template <FftScalar T = double>
class SmallDftCodelet {
public:
    explicit SmallDftCodelet(std::size_t radix) : radix_(radix) {
        if (!small_codelet_supported(radix_))
            throw std::invalid_argument("SmallDftCodelet supports radix 2, 3, 4, 5, or 7");
        if (radix_ == 3 || radix_ == 5 || radix_ == 7) {
            roots_.resize(radix_ * radix_);
            const T tau = T{2} * std::numbers::pi_v<T>;
            for (std::size_t k = 0; k < radix_; ++k)
                for (std::size_t q = 0; q < radix_; ++q)
                    roots_[k * radix_ + q] = root<T>(-tau * static_cast<T>(k * q) /
                                                      static_cast<T>(radix_));
        }
    }

    [[nodiscard]] std::size_t radix() const noexcept { return radix_; }
    [[nodiscard]] std::size_t stored_roots() const noexcept { return roots_.size(); }

    // Allocation-free, trigonometry-free execution after construction.
    void execute(std::span<ComplexT<T>> values, Direction direction) const {
        if (values.size() != radix_) throw std::invalid_argument("SmallDftCodelet buffer size mismatch");
        if (radix_ == 2) {
            const auto a = values[0], b = values[1];
            values[0] = a + b;
            values[1] = a - b;
        } else if (radix_ == 4) {
            const auto a = values[0], b = values[1], c = values[2], d = values[3];
            const auto e0 = a + c, e1 = a - c, o0 = b + d;
            const auto diff = b - d;
            const ComplexT<T> i_term = direction == Direction::Forward
                ? ComplexT<T>{diff.imag(), -diff.real()}
                : ComplexT<T>{-diff.imag(), diff.real()};
            values[0] = e0 + o0;
            values[2] = e0 - o0;
            values[1] = e1 + i_term;
            values[3] = e1 - i_term;
        } else {
            std::array<ComplexT<T>, 7> input{};
            for (std::size_t i = 0; i < radix_; ++i) input[i] = values[i];
            for (std::size_t k = 0; k < radix_; ++k) {
                ComplexT<T> sum{};
                for (std::size_t q = 0; q < radix_; ++q) {
                    auto w = roots_[k * radix_ + q];
                    if (direction == Direction::Inverse) w = std::conj(w);
                    sum += input[q] * w;
                }
                values[k] = sum;
            }
        }
        if (direction == Direction::Inverse) {
            const auto scale = T{1} / static_cast<T>(radix_);
            for (auto& z : values) z *= scale;
        }
    }

private:
    std::size_t radix_{};
    VectorT<T> roots_;
};
SmallDftCodelet(std::size_t) -> SmallDftCodelet<double>;

} // namespace fftlab
