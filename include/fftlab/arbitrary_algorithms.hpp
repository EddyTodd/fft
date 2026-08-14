#pragma once

#include "fftlab/power2_algorithms.hpp"

#include <cmath>
#include <cstddef>
#include <limits>
#include <numbers>
#include <numeric>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

namespace fftlab {

[[nodiscard]] inline std::size_t smallest_factor(std::size_t n) noexcept {
    if (n % 2 == 0) {
        return 2;
    }
    for (std::size_t factor = 3; factor <= n / factor; factor += 2) {
        if (n % factor == 0) {
            return factor;
        }
    }
    return n;
}

template <FftScalar T>
[[nodiscard]] VectorT<T> mixed_core(const VectorT<T>& input, Direction direction) {
    const auto n = input.size();
    if (n <= 1) {
        return input;
    }

    const auto radix = smallest_factor(n);
    if (radix == n) {
        auto output = dft(input, direction);
        if (direction == Direction::Inverse) {
            for (auto& value : output) {
                value *= static_cast<T>(n);
            }
        }
        return output;
    }

    const auto sub_size = n / radix;
    std::vector<VectorT<T>> subtransforms(radix, VectorT<T>(sub_size));
    for (std::size_t q = 0; q < radix; ++q) {
        for (std::size_t j = 0; j < sub_size; ++j) {
            subtransforms[q][j] = input[radix * j + q];
        }
        subtransforms[q] = mixed_core(subtransforms[q], direction);
    }

    VectorT<T> output(n);
    const T sign = sign_for<T>(direction);
    const T tau = T{2} * std::numbers::pi_v<T>;
    for (std::size_t k0 = 0; k0 < sub_size; ++k0) {
        for (std::size_t k1 = 0; k1 < radix; ++k1) {
            const auto k = k0 + sub_size * k1;
            ComplexT<T> sum{};
            for (std::size_t q = 0; q < radix; ++q) {
                const T angle = sign * tau * static_cast<T>(k) * static_cast<T>(q) /
                                static_cast<T>(n);
                sum += subtransforms[q][k0] * root<T>(angle);
            }
            output[k] = sum;
        }
    }
    return output;
}

template <FftScalar T>
[[nodiscard]] VectorT<T> mixed(const VectorT<T>& input, Direction direction) {
    auto output = mixed_core(input, direction);
    normalize_inverse<T>(output, direction);
    return output;
}

template <FftScalar T>
[[nodiscard]] VectorT<T> mixed(const VectorT<T>& input, bool inverse = false) {
    return mixed(input, inverse ? Direction::Inverse : Direction::Forward);
}

template <FftScalar T>
[[nodiscard]] VectorT<T> good_thomas(const VectorT<T>& input, Direction direction,
                                     std::size_t factor_a = 0, std::size_t factor_b = 0) {
    const auto n = input.size();
    if (n <= 1) {
        return input;
    }

    if (factor_a == 0 || factor_b == 0) {
        std::tie(factor_a, factor_b) = coprime_factor_split(n);
    }
    if (factor_a <= 1 || factor_b <= 1 || factor_a * factor_b != n ||
        std::gcd(factor_a, factor_b) != 1) {
        throw std::invalid_argument("Good-Thomas requires a coprime factorization N=a*b");
    }

    const auto inverse_b_mod_a = modular_inverse(factor_b % factor_a, factor_a);
    const auto inverse_a_mod_b = modular_inverse(factor_a % factor_b, factor_b);

    VectorT<T> matrix(n);
    for (std::size_t n1 = 0; n1 < factor_a; ++n1) {
        for (std::size_t n2 = 0; n2 < factor_b; ++n2) {
            const auto index =
                (mul_mod(mul_mod(n1, factor_b, n), inverse_b_mod_a, n) +
                 mul_mod(mul_mod(n2, factor_a, n), inverse_a_mod_b, n)) %
                n;
            matrix[n1 * factor_b + n2] = input[index];
        }
    }

    // Transform rows and columns independently; coprime indexing removes cross twiddles.
    for (std::size_t n1 = 0; n1 < factor_a; ++n1) {
        VectorT<T> row(factor_b);
        for (std::size_t n2 = 0; n2 < factor_b; ++n2) {
            row[n2] = matrix[n1 * factor_b + n2];
        }
        row = mixed(row, direction);
        for (std::size_t k2 = 0; k2 < factor_b; ++k2) {
            matrix[n1 * factor_b + k2] = row[k2];
        }
    }

    for (std::size_t k2 = 0; k2 < factor_b; ++k2) {
        VectorT<T> column(factor_a);
        for (std::size_t n1 = 0; n1 < factor_a; ++n1) {
            column[n1] = matrix[n1 * factor_b + k2];
        }
        column = mixed(column, direction);
        for (std::size_t k1 = 0; k1 < factor_a; ++k1) {
            matrix[k1 * factor_b + k2] = column[k1];
        }
    }

    VectorT<T> output(n);
    for (std::size_t k1 = 0; k1 < factor_a; ++k1) {
        for (std::size_t k2 = 0; k2 < factor_b; ++k2) {
            output[(k1 * factor_b + k2 * factor_a) % n] = matrix[k1 * factor_b + k2];
        }
    }
    // Each inverse dimension is already normalized, yielding the required 1/N total.
    return output;
}

template <FftScalar T>
[[nodiscard]] VectorT<T> good_thomas(const VectorT<T>& input, bool inverse = false,
                                     std::size_t factor_a = 0, std::size_t factor_b = 0) {
    return good_thomas(input, inverse ? Direction::Inverse : Direction::Forward, factor_a,
                       factor_b);
}

template <FftScalar T>
[[nodiscard]] VectorT<T> bluestein(const VectorT<T>& input, Direction direction) {
    const auto n = input.size();
    if (n <= 1) {
        return input;
    }
    if (n > (std::numeric_limits<std::size_t>::max() / 2) + 1) {
        throw std::length_error("Bluestein workspace overflow");
    }

    const auto convolution_size = next_pow2(2 * n - 1);
    VectorT<T> a(convolution_size);
    VectorT<T> b(convolution_size);
    const T sign = sign_for<T>(direction);

    for (std::size_t k = 0; k < n; ++k) {
        const long double kd = static_cast<long double>(k);
        const long double period = 2.0L * static_cast<long double>(n);
        const long double phase = std::fmod(kd * kd, period) / static_cast<long double>(n);
        const auto chirp = root<T>(sign * std::numbers::pi_v<T> * static_cast<T>(phase));
        a[k] = input[k] * chirp;
        b[k] = std::conj(chirp);
        if (k != 0) {
            b[convolution_size - k] = std::conj(chirp);
        }
    }

    radix2_inplace(a, Direction::Forward);
    radix2_inplace(b, Direction::Forward);
    for (std::size_t i = 0; i < convolution_size; ++i) {
        a[i] *= b[i];
    }
    radix2_inplace(a, Direction::Inverse);

    VectorT<T> output(n);
    for (std::size_t k = 0; k < n; ++k) {
        const long double kd = static_cast<long double>(k);
        const long double period = 2.0L * static_cast<long double>(n);
        const long double phase = std::fmod(kd * kd, period) / static_cast<long double>(n);
        output[k] =
            a[k] * root<T>(sign * std::numbers::pi_v<T> * static_cast<T>(phase));
    }
    normalize_inverse<T>(output, direction);
    return output;
}

template <FftScalar T>
[[nodiscard]] VectorT<T> bluestein(const VectorT<T>& input, bool inverse = false) {
    return bluestein(input, inverse ? Direction::Inverse : Direction::Forward);
}

template <FftScalar T>
[[nodiscard]] VectorT<T> circular_convolution(VectorT<T> a, VectorT<T> b) {
    if (a.size() != b.size()) {
        throw std::invalid_argument("circular convolution size mismatch");
    }
    const auto n = a.size();
    if (n == 0) {
        return {};
    }

    const auto convolution_size = next_pow2(2 * n - 1);
    a.resize(convolution_size);
    b.resize(convolution_size);
    radix2_inplace(a, Direction::Forward);
    radix2_inplace(b, Direction::Forward);
    for (std::size_t i = 0; i < convolution_size; ++i) {
        a[i] *= b[i];
    }
    radix2_inplace(a, Direction::Inverse);

    VectorT<T> output(n);
    for (std::size_t i = 0; i < 2 * n - 1; ++i) {
        output[i % n] += a[i];
    }
    return output;
}

template <FftScalar T>
[[nodiscard]] VectorT<T> rader(const VectorT<T>& input, Direction direction) {
    const auto n = input.size();
    if (n <= 2) {
        return dft(input, direction);
    }
    if (!is_prime(n)) {
        throw std::invalid_argument("Rader requires prime N");
    }

    const auto generator = primitive_root_prime(n);
    const auto cycle_size = n - 1;
    VectorT<T> a(cycle_size);
    VectorT<T> b(cycle_size);
    std::vector<std::size_t> generator_powers(cycle_size);
    generator_powers[0] = 1;
    for (std::size_t q = 1; q < cycle_size; ++q) {
        generator_powers[q] = mul_mod(generator_powers[q - 1], generator, n);
    }

    const T sign = sign_for<T>(direction);
    for (std::size_t q = 0; q < cycle_size; ++q) {
        a[q] = input[generator_powers[(cycle_size - q) % cycle_size]];
        const T angle = sign * T{2} * std::numbers::pi_v<T> *
                        static_cast<T>(generator_powers[q]) / static_cast<T>(n);
        b[q] = root<T>(angle);
    }

    const auto convolution = circular_convolution(std::move(a), std::move(b));
    VectorT<T> output(n);
    output[0] = std::accumulate(input.begin(), input.end(), ComplexT<T>{});
    for (std::size_t q = 0; q < cycle_size; ++q) {
        output[generator_powers[q]] = input[0] + convolution[q];
    }
    normalize_inverse<T>(output, direction);
    return output;
}

template <FftScalar T>
[[nodiscard]] VectorT<T> rader(const VectorT<T>& input, bool inverse = false) {
    return rader(input, inverse ? Direction::Inverse : Direction::Forward);
}

template <FftScalar T>
[[nodiscard]] VectorT<T> transform(const VectorT<T>& input,
                                   Algorithm algorithm = Algorithm::Auto,
                                   Direction direction = Direction::Forward) {
    const auto n = input.size();
    if (algorithm == Algorithm::Auto) {
        if (n <= 1) {
            return input;
        }
        if (pow2(n)) {
            algorithm = Algorithm::Radix2;
        } else if (!is_prime(n)) {
            const auto [first, second] = coprime_factor_split(n);
            algorithm = first > 1 && second > 1 ? Algorithm::GoodThomas : Algorithm::Mixed;
        } else {
            const auto rader_size = pow2(n - 1) ? n - 1 : next_pow2(2 * (n - 1) - 1);
            const auto bluestein_size = next_pow2(2 * n - 1);
            algorithm = rader_size < bluestein_size ? Algorithm::Rader : Algorithm::Bluestein;
        }
    }

    if (!supports_algorithm(algorithm, n)) {
        throw std::invalid_argument("algorithm does not support N");
    }

    switch (algorithm) {
        case Algorithm::Dft:
            return dft(input, direction);
        case Algorithm::Radix2:
            return radix2(input, direction);
        case Algorithm::Recursive:
            return recursive(input, direction);
        case Algorithm::Stockham:
            return stockham(input, direction);
        case Algorithm::Radix4:
            return radix4(input, direction);
        case Algorithm::SplitRadix:
            return split_radix(input, direction);
        case Algorithm::ModifiedSplitRadix:
            return modified_split_radix(input, direction);
        case Algorithm::Mixed:
            return mixed(input, direction);
        case Algorithm::GoodThomas:
            return good_thomas(input, direction);
        case Algorithm::Rader:
            return rader(input, direction);
        case Algorithm::Bluestein:
            return bluestein(input, direction);
        case Algorithm::Auto:
            break;
    }
    throw std::logic_error("unreachable FFT dispatch");
}

template <FftScalar T>
[[nodiscard]] VectorT<T> transform(const VectorT<T>& input, Algorithm algorithm, bool inverse) {
    return transform(input, algorithm, inverse ? Direction::Inverse : Direction::Forward);
}

}  // namespace fftlab
