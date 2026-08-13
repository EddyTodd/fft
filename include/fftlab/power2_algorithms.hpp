#pragma once

#include "fftlab/types.hpp"

namespace fftlab {

template <FftScalar T>
[[nodiscard]] inline ComplexT<T> root(T angle) {
    return {std::cos(angle), std::sin(angle)};
}

template <FftScalar T>
[[nodiscard]] inline T sign_for(Direction d) noexcept {
    return d == Direction::Forward ? T{-1} : T{1};
}

template <FftScalar T>
inline void normalize_inverse(std::span<ComplexT<T>> data, Direction d) {
    if (d != Direction::Inverse || data.empty()) return;
    const auto scale = T{1} / static_cast<T>(data.size());
    for (auto& value : data) value *= scale;
}

template <FftScalar T>
[[nodiscard]] VectorT<T> dft(const VectorT<T>& input, Direction direction) {
    const auto n = input.size();
    VectorT<T> output(n);
    if (n == 0) return output;
    const T sign = sign_for<T>(direction);
    const T tau = T{2} * std::numbers::pi_v<T>;
    for (std::size_t k = 0; k < n; ++k) {
        ComplexT<T> sum{};
        for (std::size_t t = 0; t < n; ++t) {
            const auto angle = sign * tau * static_cast<T>(k) * static_cast<T>(t) / static_cast<T>(n);
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
    if (n <= 1) return;
    if (!pow2(n)) throw std::invalid_argument("radix2 requires power-of-two N");
    for (std::size_t i = 1, j = 0; i < n; ++i) {
        std::size_t bit = n >> 1;
        for (; (j & bit) != 0; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) std::swap(data[i], data[j]);
    }
    const T sign = sign_for<T>(direction);
    const T tau = T{2} * std::numbers::pi_v<T>;
    for (std::size_t len = 2; len <= n;) {
        const auto step = root<T>(sign * tau / static_cast<T>(len));
        for (std::size_t base = 0; base < n; base += len) {
            ComplexT<T> w{1, 0};
            for (std::size_t j = 0; j < len / 2; ++j) {
                const auto u = data[base + j];
                const auto v = data[base + j + len / 2] * w;
                data[base + j] = u + v;
                data[base + j + len / 2] = u - v;
                w *= step;
            }
        }
        if (len == n) break;
        len <<= 1;
    }
    normalize_inverse<T>(data, direction);
}

template <FftScalar T>
void radix2_inplace(VectorT<T>& data, bool inverse = false) {
    radix2_inplace(data, inverse ? Direction::Inverse : Direction::Forward);
}

template <FftScalar T>
[[nodiscard]] VectorT<T> radix2(const VectorT<T>& input, Direction d) {
    auto out = input; radix2_inplace(out, d); return out;
}
template <FftScalar T>
[[nodiscard]] VectorT<T> radix2(const VectorT<T>& input, bool inverse = false) {
    return radix2(input, inverse ? Direction::Inverse : Direction::Forward);
}

template <FftScalar T>
void recursive_core(VectorT<T>& data, Direction d) {
    const auto n = data.size();
    if (n <= 1) return;
    VectorT<T> even(n / 2), odd(n / 2);
    for (std::size_t i = 0; i < n / 2; ++i) { even[i] = data[2 * i]; odd[i] = data[2 * i + 1]; }
    recursive_core(even, d); recursive_core(odd, d);
    ComplexT<T> w{1, 0};
    const auto step = root<T>(sign_for<T>(d) * T{2} * std::numbers::pi_v<T> / static_cast<T>(n));
    for (std::size_t k = 0; k < n / 2; ++k) {
        const auto v = w * odd[k];
        data[k] = even[k] + v;
        data[k + n / 2] = even[k] - v;
        w *= step;
    }
}

template <FftScalar T>
[[nodiscard]] VectorT<T> recursive(const VectorT<T>& input, Direction d) {
    if (!input.empty() && !pow2(input.size())) throw std::invalid_argument("recursive radix2 requires power-of-two N");
    auto out = input; recursive_core(out, d); normalize_inverse<T>(out, d); return out;
}
template <FftScalar T>
[[nodiscard]] VectorT<T> recursive(const VectorT<T>& input, bool inverse = false) {
    return recursive(input, inverse ? Direction::Inverse : Direction::Forward);
}

template <FftScalar T>
[[nodiscard]] VectorT<T> stockham(const VectorT<T>& input, Direction d) {
    const auto n = input.size();
    if (n <= 1) return input;
    if (!pow2(n)) throw std::invalid_argument("Stockham requires power-of-two N");
    VectorT<T> src = input, dst(n);
    const T sign = sign_for<T>(d);
    for (std::size_t m = 1; m < n; m <<= 1) {
        for (std::size_t j = 0; j < n / 2; ++j) {
            const auto k = j & (m - 1);
            const auto twiddle = root<T>(sign * std::numbers::pi_v<T> * static_cast<T>(k) / static_cast<T>(m));
            const auto v = src[j + n / 2] * twiddle;
            const auto out = 2 * j - k;
            dst[out] = src[j] + v;
            dst[out + m] = src[j] - v;
        }
        src.swap(dst);
    }
    normalize_inverse<T>(src, d);
    return src;
}
template <FftScalar T>
[[nodiscard]] VectorT<T> stockham(const VectorT<T>& input, bool inverse = false) {
    return stockham(input, inverse ? Direction::Inverse : Direction::Forward);
}

template <FftScalar T>
void radix4_core(VectorT<T>& data, Direction d) {
    const auto n = data.size();
    if (n <= 1) return;
    if (n == 2) { const auto u = data[0], v = data[1]; data[0] = u + v; data[1] = u - v; return; }
    const auto m = n / 4;
    std::array<VectorT<T>, 4> sub{VectorT<T>(m), VectorT<T>(m), VectorT<T>(m), VectorT<T>(m)};
    for (std::size_t q = 0; q < 4; ++q)
        for (std::size_t j = 0; j < m; ++j) sub[q][j] = data[4 * j + q];
    for (auto& part : sub) radix4_core(part, d);
    const T sign = sign_for<T>(d);
    const ComplexT<T> I{0, sign};
    for (std::size_t k = 0; k < m; ++k) {
        const auto b = sub[1][k] * root<T>(sign * T{2} * std::numbers::pi_v<T> * static_cast<T>(k) / static_cast<T>(n));
        const auto c = sub[2][k] * root<T>(sign * T{4} * std::numbers::pi_v<T> * static_cast<T>(k) / static_cast<T>(n));
        const auto e = sub[3][k] * root<T>(sign * T{6} * std::numbers::pi_v<T> * static_cast<T>(k) / static_cast<T>(n));
        const auto t0 = sub[0][k] + c, t1 = sub[0][k] - c;
        const auto t2 = b + e, t3 = (b - e) * I;
        data[k] = t0 + t2; data[k + m] = t1 + t3;
        data[k + 2 * m] = t0 - t2; data[k + 3 * m] = t1 - t3;
    }
}

template <FftScalar T>
[[nodiscard]] VectorT<T> radix4(const VectorT<T>& input, Direction d) {
    if (!input.empty() && !pow2(input.size())) throw std::invalid_argument("radix4 requires power-of-two N");
    auto out = input; radix4_core(out, d); normalize_inverse<T>(out, d); return out;
}
template <FftScalar T>
[[nodiscard]] VectorT<T> radix4(const VectorT<T>& input, bool inverse = false) {
    return radix4(input, inverse ? Direction::Inverse : Direction::Forward);
}

template <FftScalar T>
void split_core(VectorT<T>& data, Direction d) {
    const auto n = data.size();
    if (n <= 1) return;
    if (n == 2) { const auto u = data[0], v = data[1]; data[0] = u + v; data[1] = u - v; return; }
    const auto qn = n / 4;
    VectorT<T> even(n / 2), odd1(qn), odd3(qn);
    for (std::size_t j = 0; j < n / 2; ++j) even[j] = data[2 * j];
    for (std::size_t j = 0; j < qn; ++j) { odd1[j] = data[4 * j + 1]; odd3[j] = data[4 * j + 3]; }
    split_core(even, d); split_core(odd1, d); split_core(odd3, d);
    const T sign = sign_for<T>(d);
    const ComplexT<T> I{0, sign};
    for (std::size_t k = 0; k < qn; ++k) {
        const auto t1 = odd1[k] * root<T>(sign * T{2} * std::numbers::pi_v<T> * static_cast<T>(k) / static_cast<T>(n));
        const auto t2 = odd3[k] * root<T>(sign * T{6} * std::numbers::pi_v<T> * static_cast<T>(k) / static_cast<T>(n));
        const auto s = t1 + t2, diff = t1 - t2;
        data[k] = even[k] + s; data[k + n / 2] = even[k] - s;
        data[k + qn] = even[k + qn] + I * diff;
        data[k + 3 * qn] = even[k + qn] - I * diff;
    }
}

template <FftScalar T>
[[nodiscard]] VectorT<T> split_radix(const VectorT<T>& input, Direction d) {
    if (!input.empty() && !pow2(input.size())) throw std::invalid_argument("split-radix requires power-of-two N");
    auto out = input; split_core(out, d); normalize_inverse<T>(out, d); return out;
}
template <FftScalar T>
[[nodiscard]] VectorT<T> split_radix(const VectorT<T>& input, bool inverse = false) {
    return split_radix(input, inverse ? Direction::Inverse : Direction::Forward);
}

// Johnson-Frigo scale schedule. The implementation below is a correctness-first
// scaled conjugate-pair split-radix treatment: it uses the published recursive
// scale family and its 2-real-multiply modified twiddle forms. It is intentionally
// readable rather than a generated flop-minimal codelet implementation.
template <FftScalar T>
[[nodiscard]] T modified_split_scale(std::size_t n, std::size_t k) {
    if (n <= 4) return T{1};
    const auto quarter = n / 4;
    const auto kk = k % quarter;
    const auto sub = modified_split_scale<T>(quarter, kk);
    const auto angle = T{2} * std::numbers::pi_v<T> * static_cast<T>(kk) / static_cast<T>(n);
    return sub * (kk <= n / 8 ? std::cos(angle) : std::sin(angle));
}

template <FftScalar T>
[[nodiscard]] ComplexT<T> modified_twiddle_mul(const ComplexT<T>& z, std::size_t n,
                                               std::size_t k, Direction d, bool conjugate) {
    const auto kk = k % (n / 4);
    const auto theta = T{2} * std::numbers::pi_v<T> * static_cast<T>(kk) / static_cast<T>(n);
    // forward t is either 1-i*tan(theta) or cot(theta)-i. Inverse/conjugate flip signs.
    T imag_sign = d == Direction::Forward ? T{-1} : T{1};
    if (conjugate) imag_sign = -imag_sign;
    const T a = z.real(), b = z.imag();
    if (kk <= n / 8) {
        const T t = std::tan(theta);
        return {a - imag_sign * b * t, b + imag_sign * a * t};
    }
    const T cot = T{1} / std::tan(theta);
    return {a * cot - imag_sign * b, b * cot + imag_sign * a};
}

template <FftScalar T>
void modified_split_scaled_core(const VectorT<T>& input, VectorT<T>& scaled, Direction d) {
    const auto n = input.size();
    scaled.resize(n);
    if (n <= 4) {
        scaled = dft(input, d);
        // dft() normalizes inverse; recursive modified split must remain unnormalized.
        if (d == Direction::Inverse && n > 0) {
            for (auto& z : scaled) z *= static_cast<T>(n);
        }
        return;
    }
    VectorT<T> even_input(n / 2), odd1_input(n / 4), odd3_input(n / 4);
    for (std::size_t j = 0; j < n / 2; ++j) even_input[j] = input[2 * j];
    for (std::size_t j = 0; j < n / 4; ++j) {
        odd1_input[j] = input[4 * j + 1];
        odd3_input[j] = input[(4 * j + n - 1) % n];
    }
    VectorT<T> even, odd1, odd3;
    modified_split_scaled_core(even_input, even, d);
    modified_split_scaled_core(odd1_input, odd1, d);
    modified_split_scaled_core(odd3_input, odd3, d);
    const auto quarter = n / 4;
    for (std::size_t k = 0; k < quarter; ++k) {
        const T scale_n = modified_split_scale<T>(n, k);
        const auto t1 = modified_twiddle_mul<T>(odd1[k], n, k, d, false);
        const auto t2 = modified_twiddle_mul<T>(odd3[k], n, k, d, true);
        const auto sum = t1 + t2, diff = t1 - t2;
        const auto e0 = even[k] * (modified_split_scale<T>(n / 2, k) / scale_n);
        const auto e1 = even[k + quarter] * (modified_split_scale<T>(n / 2, k + quarter) / scale_n);
        const ComplexT<T> minus_i = d == Direction::Forward ? ComplexT<T>{0, -1} : ComplexT<T>{0, 1};
        scaled[k] = e0 + sum;
        scaled[k + quarter] = e1 + minus_i * diff;
        scaled[k + 2 * quarter] = e0 - sum;
        scaled[k + 3 * quarter] = e1 - minus_i * diff;
    }
}

template <FftScalar T>
[[nodiscard]] VectorT<T> modified_split_radix(const VectorT<T>& input, Direction d) {
    if (!input.empty() && !pow2(input.size())) throw std::invalid_argument("modified split-radix requires power-of-two N");
    if (input.size() <= 1) return input;
    VectorT<T> scaled;
    modified_split_scaled_core(input, scaled, d);
    for (std::size_t k = 0; k < scaled.size(); ++k)
        scaled[k] *= modified_split_scale<T>(scaled.size(), k);
    normalize_inverse<T>(scaled, d);
    return scaled;
}
template <FftScalar T>
[[nodiscard]] VectorT<T> modified_split_radix(const VectorT<T>& input, bool inverse = false) {
    return modified_split_radix(input, inverse ? Direction::Inverse : Direction::Forward);
}


} // namespace fftlab
