#pragma once

#include "fftlab/types.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <span>
#include <stdexcept>
#include <utility>

namespace fftlab {

template <FftScalar T>
[[nodiscard]] inline ComplexT<T> root(T angle) {
    return {std::cos(angle), std::sin(angle)};
}

template <FftScalar T>
[[nodiscard]] inline T sign_for(Direction direction) noexcept {
    return direction == Direction::Forward ? T{-1} : T{1};
}

template <FftScalar T>
inline void normalize_inverse(std::span<ComplexT<T>> data, Direction direction) {
    if (direction != Direction::Inverse || data.empty()) {
        return;
    }
    const auto scale = T{1} / static_cast<T>(data.size());
    for (auto& value : data) {
        value *= scale;
    }
}

template <FftScalar T>
[[nodiscard]] VectorT<T> dft(const VectorT<T>& input, Direction direction) {
    const auto n = input.size();
    VectorT<T> output(n);
    if (n == 0) {
        return output;
    }

    const T sign = sign_for<T>(direction);
    const T tau = T{2} * std::numbers::pi_v<T>;
    for (std::size_t k = 0; k < n; ++k) {
        ComplexT<T> sum{};
        for (std::size_t t = 0; t < n; ++t) {
            const auto angle = sign * tau * static_cast<T>(k) * static_cast<T>(t) /
                               static_cast<T>(n);
            sum += input[t] * root<T>(angle);
        }
        output[k] = sum;
    }
    normalize_inverse<T>(output, direction);
    return output;
}

template <FftScalar T>
[[nodiscard]] VectorT<T> dft(const VectorT<T>& input, bool inverse = false) {
    return dft(input, inverse ? Direction::Inverse : Direction::Forward);
}

template <FftScalar T>
void radix2_inplace(VectorT<T>& data, Direction direction) {
    const auto n = data.size();
    if (n <= 1) {
        return;
    }
    if (!pow2(n)) {
        throw std::invalid_argument("radix2 requires power-of-two N");
    }

    for (std::size_t i = 1, j = 0; i < n; ++i) {
        std::size_t bit = n >> 1;
        for (; (j & bit) != 0; bit >>= 1) {
            j ^= bit;
        }
        j ^= bit;
        if (i < j) {
            std::swap(data[i], data[j]);
        }
    }

    const T sign = sign_for<T>(direction);
    const T tau = T{2} * std::numbers::pi_v<T>;
    for (std::size_t length = 2; length <= n;) {
        const auto step = root<T>(sign * tau / static_cast<T>(length));
        for (std::size_t base = 0; base < n; base += length) {
            ComplexT<T> twiddle{1, 0};
            for (std::size_t offset = 0; offset < length / 2; ++offset) {
                const auto upper = data[base + offset];
                const auto lower = data[base + offset + length / 2] * twiddle;
                data[base + offset] = upper + lower;
                data[base + offset + length / 2] = upper - lower;
                twiddle *= step;
            }
        }
        if (length == n) {
            break;
        }
        length <<= 1;
    }
    normalize_inverse<T>(data, direction);
}

template <FftScalar T>
void radix2_inplace(VectorT<T>& data, bool inverse = false) {
    radix2_inplace(data, inverse ? Direction::Inverse : Direction::Forward);
}

template <FftScalar T>
[[nodiscard]] VectorT<T> radix2(const VectorT<T>& input, Direction direction) {
    auto output = input;
    radix2_inplace(output, direction);
    return output;
}

template <FftScalar T>
[[nodiscard]] VectorT<T> radix2(const VectorT<T>& input, bool inverse = false) {
    return radix2(input, inverse ? Direction::Inverse : Direction::Forward);
}

template <FftScalar T>
void recursive_core(VectorT<T>& data, Direction direction) {
    const auto n = data.size();
    if (n <= 1) {
        return;
    }

    VectorT<T> even(n / 2);
    VectorT<T> odd(n / 2);
    for (std::size_t i = 0; i < n / 2; ++i) {
        even[i] = data[2 * i];
        odd[i] = data[2 * i + 1];
    }

    recursive_core(even, direction);
    recursive_core(odd, direction);

    ComplexT<T> twiddle{1, 0};
    const auto step =
        root<T>(sign_for<T>(direction) * T{2} * std::numbers::pi_v<T> / static_cast<T>(n));
    for (std::size_t k = 0; k < n / 2; ++k) {
        const auto lower = twiddle * odd[k];
        data[k] = even[k] + lower;
        data[k + n / 2] = even[k] - lower;
        twiddle *= step;
    }
}

template <FftScalar T>
[[nodiscard]] VectorT<T> recursive(const VectorT<T>& input, Direction direction) {
    if (!input.empty() && !pow2(input.size())) {
        throw std::invalid_argument("recursive radix2 requires power-of-two N");
    }
    auto output = input;
    recursive_core(output, direction);
    normalize_inverse<T>(output, direction);
    return output;
}

template <FftScalar T>
[[nodiscard]] VectorT<T> recursive(const VectorT<T>& input, bool inverse = false) {
    return recursive(input, inverse ? Direction::Inverse : Direction::Forward);
}

template <FftScalar T>
[[nodiscard]] VectorT<T> stockham(const VectorT<T>& input, Direction direction) {
    const auto n = input.size();
    if (n <= 1) {
        return input;
    }
    if (!pow2(n)) {
        throw std::invalid_argument("Stockham requires power-of-two N");
    }

    VectorT<T> source = input;
    VectorT<T> destination(n);
    const T sign = sign_for<T>(direction);
    for (std::size_t m = 1; m < n; m <<= 1) {
        for (std::size_t j = 0; j < n / 2; ++j) {
            const auto k = j & (m - 1);
            const auto twiddle = root<T>(sign * std::numbers::pi_v<T> * static_cast<T>(k) /
                                         static_cast<T>(m));
            const auto lower = source[j + n / 2] * twiddle;
            const auto output_index = 2 * j - k;
            destination[output_index] = source[j] + lower;
            destination[output_index + m] = source[j] - lower;
        }
        source.swap(destination);
    }
    normalize_inverse<T>(source, direction);
    return source;
}

template <FftScalar T>
[[nodiscard]] VectorT<T> stockham(const VectorT<T>& input, bool inverse = false) {
    return stockham(input, inverse ? Direction::Inverse : Direction::Forward);
}

template <FftScalar T>
void radix4_core(VectorT<T>& data, Direction direction) {
    const auto n = data.size();
    if (n <= 1) {
        return;
    }
    if (n == 2) {
        const auto first = data[0];
        const auto second = data[1];
        data[0] = first + second;
        data[1] = first - second;
        return;
    }

    const auto quarter = n / 4;
    std::array<VectorT<T>, 4> subtransforms{
        VectorT<T>(quarter),
        VectorT<T>(quarter),
        VectorT<T>(quarter),
        VectorT<T>(quarter),
    };
    for (std::size_t residue = 0; residue < 4; ++residue) {
        for (std::size_t j = 0; j < quarter; ++j) {
            subtransforms[residue][j] = data[4 * j + residue];
        }
    }
    for (auto& subtransform : subtransforms) {
        radix4_core(subtransform, direction);
    }

    const T sign = sign_for<T>(direction);
    const ComplexT<T> imaginary_unit{0, sign};
    const T tau = T{2} * std::numbers::pi_v<T>;
    for (std::size_t k = 0; k < quarter; ++k) {
        const auto b = subtransforms[1][k] *
                       root<T>(sign * tau * static_cast<T>(k) / static_cast<T>(n));
        const auto c = subtransforms[2][k] *
                       root<T>(sign * T{2} * tau * static_cast<T>(k) / static_cast<T>(n));
        const auto d = subtransforms[3][k] *
                       root<T>(sign * T{3} * tau * static_cast<T>(k) / static_cast<T>(n));
        const auto sum_even = subtransforms[0][k] + c;
        const auto difference_even = subtransforms[0][k] - c;
        const auto sum_odd = b + d;
        const auto difference_odd = (b - d) * imaginary_unit;

        data[k] = sum_even + sum_odd;
        data[k + quarter] = difference_even + difference_odd;
        data[k + 2 * quarter] = sum_even - sum_odd;
        data[k + 3 * quarter] = difference_even - difference_odd;
    }
}

template <FftScalar T>
[[nodiscard]] VectorT<T> radix4(const VectorT<T>& input, Direction direction) {
    if (!input.empty() && !pow2(input.size())) {
        throw std::invalid_argument("radix4 requires power-of-two N");
    }
    auto output = input;
    radix4_core(output, direction);
    normalize_inverse<T>(output, direction);
    return output;
}

template <FftScalar T>
[[nodiscard]] VectorT<T> radix4(const VectorT<T>& input, bool inverse = false) {
    return radix4(input, inverse ? Direction::Inverse : Direction::Forward);
}

template <FftScalar T>
void split_core(VectorT<T>& data, Direction direction) {
    const auto n = data.size();
    if (n <= 1) {
        return;
    }
    if (n == 2) {
        const auto first = data[0];
        const auto second = data[1];
        data[0] = first + second;
        data[1] = first - second;
        return;
    }

    const auto quarter = n / 4;
    VectorT<T> even(n / 2);
    VectorT<T> odd_one(quarter);
    VectorT<T> odd_three(quarter);
    for (std::size_t j = 0; j < n / 2; ++j) {
        even[j] = data[2 * j];
    }
    for (std::size_t j = 0; j < quarter; ++j) {
        odd_one[j] = data[4 * j + 1];
        odd_three[j] = data[4 * j + 3];
    }

    split_core(even, direction);
    split_core(odd_one, direction);
    split_core(odd_three, direction);

    const T sign = sign_for<T>(direction);
    const ComplexT<T> imaginary_unit{0, sign};
    const T tau = T{2} * std::numbers::pi_v<T>;
    for (std::size_t k = 0; k < quarter; ++k) {
        const auto first_twiddle =
            odd_one[k] * root<T>(sign * tau * static_cast<T>(k) / static_cast<T>(n));
        const auto third_twiddle = odd_three[k] *
                                   root<T>(sign * T{3} * tau * static_cast<T>(k) /
                                           static_cast<T>(n));
        const auto sum = first_twiddle + third_twiddle;
        const auto difference = first_twiddle - third_twiddle;

        data[k] = even[k] + sum;
        data[k + n / 2] = even[k] - sum;
        data[k + quarter] = even[k + quarter] + imaginary_unit * difference;
        data[k + 3 * quarter] = even[k + quarter] - imaginary_unit * difference;
    }
}

template <FftScalar T>
[[nodiscard]] VectorT<T> split_radix(const VectorT<T>& input, Direction direction) {
    if (!input.empty() && !pow2(input.size())) {
        throw std::invalid_argument("split-radix requires power-of-two N");
    }
    auto output = input;
    split_core(output, direction);
    normalize_inverse<T>(output, direction);
    return output;
}

template <FftScalar T>
[[nodiscard]] VectorT<T> split_radix(const VectorT<T>& input, bool inverse = false) {
    return split_radix(input, inverse ? Direction::Inverse : Direction::Forward);
}

// Johnson-Frigo scale schedule. This correctness-first implementation uses the
// published recursive scale family and two-real-multiply modified twiddle forms.
template <FftScalar T>
[[nodiscard]] T modified_split_scale(std::size_t n, std::size_t k) {
    if (n <= 4) {
        return T{1};
    }
    const auto quarter = n / 4;
    const auto reduced_k = k % quarter;
    const auto subscale = modified_split_scale<T>(quarter, reduced_k);
    const auto angle = T{2} * std::numbers::pi_v<T> * static_cast<T>(reduced_k) /
                       static_cast<T>(n);
    return subscale * (reduced_k <= n / 8 ? std::cos(angle) : std::sin(angle));
}

template <FftScalar T>
[[nodiscard]] ComplexT<T> modified_twiddle_mul(const ComplexT<T>& value, std::size_t n,
                                               std::size_t k, Direction direction,
                                               bool conjugate) {
    const auto reduced_k = k % (n / 4);
    const auto theta = T{2} * std::numbers::pi_v<T> * static_cast<T>(reduced_k) /
                       static_cast<T>(n);

    T imaginary_sign = direction == Direction::Forward ? T{-1} : T{1};
    if (conjugate) {
        imaginary_sign = -imaginary_sign;
    }

    const T real = value.real();
    const T imaginary = value.imag();
    if (reduced_k <= n / 8) {
        const T tangent = std::tan(theta);
        return {real - imaginary_sign * imaginary * tangent,
                imaginary + imaginary_sign * real * tangent};
    }

    const T cotangent = T{1} / std::tan(theta);
    return {real * cotangent - imaginary_sign * imaginary,
            imaginary * cotangent + imaginary_sign * real};
}

template <FftScalar T>
void modified_split_scaled_core(const VectorT<T>& input, VectorT<T>& scaled,
                                Direction direction) {
    const auto n = input.size();
    scaled.resize(n);
    if (n <= 4) {
        scaled = dft(input, direction);
        // dft() normalizes inverse transforms; recursive split-radix must remain unnormalized.
        if (direction == Direction::Inverse && n > 0) {
            for (auto& value : scaled) {
                value *= static_cast<T>(n);
            }
        }
        return;
    }

    VectorT<T> even_input(n / 2);
    VectorT<T> odd_one_input(n / 4);
    VectorT<T> odd_three_input(n / 4);
    for (std::size_t j = 0; j < n / 2; ++j) {
        even_input[j] = input[2 * j];
    }
    for (std::size_t j = 0; j < n / 4; ++j) {
        odd_one_input[j] = input[4 * j + 1];
        odd_three_input[j] = input[(4 * j + n - 1) % n];
    }

    VectorT<T> even;
    VectorT<T> odd_one;
    VectorT<T> odd_three;
    modified_split_scaled_core(even_input, even, direction);
    modified_split_scaled_core(odd_one_input, odd_one, direction);
    modified_split_scaled_core(odd_three_input, odd_three, direction);

    const auto quarter = n / 4;
    for (std::size_t k = 0; k < quarter; ++k) {
        const T scale = modified_split_scale<T>(n, k);
        const auto first_twiddle = modified_twiddle_mul<T>(odd_one[k], n, k, direction, false);
        const auto third_twiddle =
            modified_twiddle_mul<T>(odd_three[k], n, k, direction, true);
        const auto sum = first_twiddle + third_twiddle;
        const auto difference = first_twiddle - third_twiddle;
        const auto even_zero = even[k] * (modified_split_scale<T>(n / 2, k) / scale);
        const auto even_one =
            even[k + quarter] * (modified_split_scale<T>(n / 2, k + quarter) / scale);
        const ComplexT<T> minus_i =
            direction == Direction::Forward ? ComplexT<T>{0, -1} : ComplexT<T>{0, 1};

        scaled[k] = even_zero + sum;
        scaled[k + quarter] = even_one + minus_i * difference;
        scaled[k + 2 * quarter] = even_zero - sum;
        scaled[k + 3 * quarter] = even_one - minus_i * difference;
    }
}

template <FftScalar T>
[[nodiscard]] VectorT<T> modified_split_radix(const VectorT<T>& input, Direction direction) {
    if (!input.empty() && !pow2(input.size())) {
        throw std::invalid_argument("modified split-radix requires power-of-two N");
    }
    if (input.size() <= 1) {
        return input;
    }

    VectorT<T> scaled;
    modified_split_scaled_core(input, scaled, direction);
    for (std::size_t k = 0; k < scaled.size(); ++k) {
        scaled[k] *= modified_split_scale<T>(scaled.size(), k);
    }
    normalize_inverse<T>(scaled, direction);
    return scaled;
}

template <FftScalar T>
[[nodiscard]] VectorT<T> modified_split_radix(const VectorT<T>& input, bool inverse = false) {
    return modified_split_radix(input, inverse ? Direction::Inverse : Direction::Forward);
}

}  // namespace fftlab
