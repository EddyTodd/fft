#pragma once

#include "fftlab/power2_algorithms.hpp"

namespace fftlab {

[[nodiscard]] inline std::size_t smallest_factor(std::size_t n) noexcept {
    if (n % 2 == 0) return 2;
    for (std::size_t p = 3; p <= n / p; p += 2) if (n % p == 0) return p;
    return n;
}

template <FftScalar T>
[[nodiscard]] VectorT<T> mixed_core(const VectorT<T>& input, Direction d) {
    const auto n = input.size();
    if (n <= 1) return input;
    const auto r = smallest_factor(n);
    if (r == n) {
        auto out = dft(input, d);
        if (d == Direction::Inverse) for (auto& z : out) z *= static_cast<T>(n);
        return out;
    }
    const auto m = n / r;
    std::vector<VectorT<T>> sub(r, VectorT<T>(m));
    for (std::size_t q = 0; q < r; ++q) {
        for (std::size_t j = 0; j < m; ++j) sub[q][j] = input[r * j + q];
        sub[q] = mixed_core(sub[q], d);
    }
    VectorT<T> output(n);
    const T sign = sign_for<T>(d);
    const T tau = T{2} * std::numbers::pi_v<T>;
    for (std::size_t k0 = 0; k0 < m; ++k0) {
        for (std::size_t k1 = 0; k1 < r; ++k1) {
            const auto k = k0 + m * k1;
            ComplexT<T> sum{};
            for (std::size_t q = 0; q < r; ++q) {
                sum += sub[q][k0] * root<T>(sign * tau * static_cast<T>(k) * static_cast<T>(q) / static_cast<T>(n));
            }
            output[k] = sum;
        }
    }
    return output;
}

template <FftScalar T>
[[nodiscard]] VectorT<T> mixed(const VectorT<T>& input, Direction d) {
    auto out = mixed_core(input, d);
    normalize_inverse<T>(out, d);
    return out;
}
template <FftScalar T>
[[nodiscard]] VectorT<T> mixed(const VectorT<T>& input, bool inverse = false) {
    return mixed(input, inverse ? Direction::Inverse : Direction::Forward);
}

template <FftScalar T>
[[nodiscard]] VectorT<T> good_thomas(const VectorT<T>& input, Direction d,
                                     std::size_t factor_a = 0, std::size_t factor_b = 0) {
    const auto n = input.size();
    if (n <= 1) return input;
    if (factor_a == 0 || factor_b == 0) {
        std::tie(factor_a, factor_b) = coprime_factor_split(n);
    }
    if (factor_a <= 1 || factor_b <= 1 || factor_a * factor_b != n || std::gcd(factor_a, factor_b) != 1)
        throw std::invalid_argument("Good-Thomas requires a coprime factorization N=a*b");
    const auto inv_b_mod_a = modular_inverse(factor_b % factor_a, factor_a);
    const auto inv_a_mod_b = modular_inverse(factor_a % factor_b, factor_b);
    VectorT<T> matrix(n);
    for (std::size_t n1 = 0; n1 < factor_a; ++n1) {
        for (std::size_t n2 = 0; n2 < factor_b; ++n2) {
            const auto index = (mul_mod(mul_mod(n1, factor_b, n), inv_b_mod_a, n) +
                                mul_mod(mul_mod(n2, factor_a, n), inv_a_mod_b, n)) % n;
            matrix[n1 * factor_b + n2] = input[index];
        }
    }
    // Transform rows (factor_b) and columns (factor_a); no cross twiddles appear.
    for (std::size_t n1 = 0; n1 < factor_a; ++n1) {
        VectorT<T> row(factor_b);
        for (std::size_t n2 = 0; n2 < factor_b; ++n2) row[n2] = matrix[n1 * factor_b + n2];
        row = mixed(row, d);
        for (std::size_t k2 = 0; k2 < factor_b; ++k2) matrix[n1 * factor_b + k2] = row[k2];
    }
    for (std::size_t k2 = 0; k2 < factor_b; ++k2) {
        VectorT<T> col(factor_a);
        for (std::size_t n1 = 0; n1 < factor_a; ++n1) col[n1] = matrix[n1 * factor_b + k2];
        col = mixed(col, d);
        for (std::size_t k1 = 0; k1 < factor_a; ++k1) matrix[k1 * factor_b + k2] = col[k1];
    }
    VectorT<T> output(n);
    for (std::size_t k1 = 0; k1 < factor_a; ++k1)
        for (std::size_t k2 = 0; k2 < factor_b; ++k2)
            output[(k1 * factor_b + k2 * factor_a) % n] = matrix[k1 * factor_b + k2];
    // Each dimension's inverse was already normalized, yielding 1/N total.
    return output;
}
template <FftScalar T>
[[nodiscard]] VectorT<T> good_thomas(const VectorT<T>& input, bool inverse = false,
                                     std::size_t factor_a = 0, std::size_t factor_b = 0) {
    return good_thomas(input, inverse ? Direction::Inverse : Direction::Forward, factor_a, factor_b);
}

template <FftScalar T>
[[nodiscard]] VectorT<T> bluestein(const VectorT<T>& input, Direction d) {
    const auto n = input.size();
    if (n <= 1) return input;
    if (n > (std::numeric_limits<std::size_t>::max() / 2) + 1)
        throw std::length_error("Bluestein workspace overflow");
    const auto m = next_pow2(2 * n - 1);
    VectorT<T> a(m), b(m);
    const T sign = sign_for<T>(d);
    for (std::size_t k = 0; k < n; ++k) {
        const long double kd = static_cast<long double>(k);
        const long double period = 2.0L * static_cast<long double>(n);
        const long double phase = std::fmod(kd * kd, period) / static_cast<long double>(n);
        const auto c = root<T>(sign * std::numbers::pi_v<T> * static_cast<T>(phase));
        a[k] = input[k] * c;
        b[k] = std::conj(c);
        if (k != 0) b[m - k] = std::conj(c);
    }
    radix2_inplace(a, Direction::Forward);
    radix2_inplace(b, Direction::Forward);
    for (std::size_t i = 0; i < m; ++i) a[i] *= b[i];
    radix2_inplace(a, Direction::Inverse);
    VectorT<T> output(n);
    for (std::size_t k = 0; k < n; ++k) {
        const long double kd = static_cast<long double>(k);
        const long double period = 2.0L * static_cast<long double>(n);
        const long double phase = std::fmod(kd * kd, period) / static_cast<long double>(n);
        output[k] = a[k] * root<T>(sign * std::numbers::pi_v<T> * static_cast<T>(phase));
    }
    normalize_inverse<T>(output, d);
    return output;
}
template <FftScalar T>
[[nodiscard]] VectorT<T> bluestein(const VectorT<T>& input, bool inverse = false) {
    return bluestein(input, inverse ? Direction::Inverse : Direction::Forward);
}

template <FftScalar T>
[[nodiscard]] VectorT<T> circular_convolution(VectorT<T> a, VectorT<T> b) {
    if (a.size() != b.size()) throw std::invalid_argument("circular convolution size mismatch");
    const auto n = a.size();
    if (n == 0) return {};
    const auto m = next_pow2(2 * n - 1);
    a.resize(m); b.resize(m);
    radix2_inplace(a, Direction::Forward); radix2_inplace(b, Direction::Forward);
    for (std::size_t i = 0; i < m; ++i) a[i] *= b[i];
    radix2_inplace(a, Direction::Inverse);
    VectorT<T> out(n);
    for (std::size_t i = 0; i < 2 * n - 1; ++i) out[i % n] += a[i];
    return out;
}

template <FftScalar T>
[[nodiscard]] VectorT<T> rader(const VectorT<T>& input, Direction d) {
    const auto n = input.size();
    if (n <= 2) return dft(input, d);
    if (!is_prime(n)) throw std::invalid_argument("Rader requires prime N");
    const auto g = primitive_root_prime(n), l = n - 1;
    VectorT<T> a(l), b(l);
    std::vector<std::size_t> gp(l);
    gp[0] = 1;
    for (std::size_t q = 1; q < l; ++q) gp[q] = mul_mod(gp[q - 1], g, n);
    const T sign = sign_for<T>(d);
    for (std::size_t q = 0; q < l; ++q) {
        a[q] = input[gp[(l - q) % l]];
        b[q] = root<T>(sign * T{2} * std::numbers::pi_v<T> * static_cast<T>(gp[q]) / static_cast<T>(n));
    }
    const auto c = circular_convolution(std::move(a), std::move(b));
    VectorT<T> output(n);
    output[0] = std::accumulate(input.begin(), input.end(), ComplexT<T>{});
    for (std::size_t q = 0; q < l; ++q) output[gp[q]] = input[0] + c[q];
    normalize_inverse<T>(output, d);
    return output;
}
template <FftScalar T>
[[nodiscard]] VectorT<T> rader(const VectorT<T>& input, bool inverse = false) {
    return rader(input, inverse ? Direction::Inverse : Direction::Forward);
}

template <FftScalar T>
[[nodiscard]] VectorT<T> transform(const VectorT<T>& input, Algo a = Algo::Auto,
                                   Direction d = Direction::Forward) {
    const auto n = input.size();
    if (a == Algo::Auto) {
        if (n <= 1) return input;
        if (pow2(n)) a = Algo::Radix2;
        else if (!is_prime(n)) {
            const auto [x, y] = coprime_factor_split(n);
            a = (x > 1 && y > 1) ? Algo::GoodThomas : Algo::Mixed;
        } else {
            const auto rader_m = pow2(n - 1) ? n - 1 : next_pow2(2 * (n - 1) - 1);
            const auto blue_m = next_pow2(2 * n - 1);
            a = rader_m < blue_m ? Algo::Rader : Algo::Bluestein;
        }
    }
    if (!supports(a, n)) throw std::invalid_argument("algorithm does not support N");
    switch (a) {
        case Algo::Dft: return dft(input, d);
        case Algo::Radix2: return radix2(input, d);
        case Algo::Recursive: return recursive(input, d);
        case Algo::Stockham: return stockham(input, d);
        case Algo::Radix4: return radix4(input, d);
        case Algo::SplitRadix: return split_radix(input, d);
        case Algo::ModifiedSplitRadix: return modified_split_radix(input, d);
        case Algo::Mixed: return mixed(input, d);
        case Algo::GoodThomas: return good_thomas(input, d);
        case Algo::Rader: return rader(input, d);
        case Algo::Bluestein: return bluestein(input, d);
        case Algo::Auto: break;
    }
    throw std::logic_error("unreachable FFT dispatch");
}

template <FftScalar T>
[[nodiscard]] VectorT<T> transform(const VectorT<T>& input, Algo a, bool inverse) {
    return transform(input, a, inverse ? Direction::Inverse : Direction::Forward);
}


} // namespace fftlab
