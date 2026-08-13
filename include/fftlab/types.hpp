#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numbers>
#include <numeric>
#include <span>
#include <stdexcept>
#include <string_view>
#include <tuple>
#include <utility>
#include <vector>

namespace fftlab {

template <class T>
concept FftScalar = std::same_as<T, float> || std::same_as<T, double>;

template <FftScalar T> using ComplexT = std::complex<T>;
template <FftScalar T> using VectorT = std::vector<ComplexT<T>>;
using Complex32 = ComplexT<float>;
using Complex64 = ComplexT<double>;
using Vector32 = VectorT<float>;
using Vector64 = VectorT<double>;
// Compatibility aliases for the historical binary64 API.
using Complex = Complex64;
using Vector = Vector64;

enum class Direction { Forward, Inverse };
enum class Algo {
    Auto,
    Dft,
    Radix2,
    Recursive,
    Stockham,
    Radix4,
    SplitRadix,
    ModifiedSplitRadix,
    Mixed,
    GoodThomas,
    Rader,
    Bluestein
};
using Algorithm = Algo;
inline constexpr std::array<Algo, 12> all_algos{
    Algo::Auto, Algo::Dft, Algo::Radix2, Algo::Recursive, Algo::Stockham, Algo::Radix4,
    Algo::SplitRadix, Algo::ModifiedSplitRadix, Algo::Mixed, Algo::GoodThomas,
    Algo::Rader, Algo::Bluestein};

[[nodiscard]] inline constexpr bool pow2(std::size_t n) noexcept {
    return n != 0 && (n & (n - 1)) == 0;
}

[[nodiscard]] inline bool is_prime(std::size_t n) noexcept {
    if (n < 2) return false;
    if ((n & 1U) == 0U) return n == 2;
    for (std::size_t d = 3; d <= n / d; d += 2) {
        if (n % d == 0) return false;
    }
    return true;
}

[[nodiscard]] inline std::size_t next_pow2(std::size_t n) {
    if (n <= 1) return 1;
    const auto high = std::size_t{1} << (std::numeric_limits<std::size_t>::digits - 1);
    if (n > high) throw std::overflow_error("next power of two overflows size_t");
    --n;
    for (std::size_t shift = 1; shift < std::numeric_limits<std::size_t>::digits; shift <<= 1)
        n |= n >> shift;
    return n + 1;
}

[[nodiscard]] inline std::size_t ilog2(std::size_t n) {
    if (!pow2(n)) throw std::invalid_argument("ilog2 requires a power of two");
    std::size_t out = 0;
    while (n > 1) { n >>= 1; ++out; }
    return out;
}

[[nodiscard]] inline std::size_t mul_mod(std::size_t a, std::size_t b, std::size_t m) noexcept {
    std::size_t result = 0;
    while (b != 0) {
        if (b & 1U) result = a >= m - result ? a - (m - result) : a + result;
        b >>= 1;
        if (b != 0) a = a >= m - a ? a - (m - a) : a + a;
    }
    return result;
}

[[nodiscard]] inline std::size_t pow_mod(std::size_t a, std::size_t e, std::size_t m) noexcept {
    std::size_t result = 1 % m;
    while (e != 0) {
        if (e & 1U) result = mul_mod(result, a, m);
        e >>= 1;
        if (e != 0) a = mul_mod(a, a, m);
    }
    return result;
}

[[nodiscard]] inline std::vector<std::size_t> prime_factors(std::size_t n) {
    std::vector<std::size_t> out;
    for (std::size_t p = 2; p <= n / p; ++p) {
        if (n % p != 0) continue;
        out.push_back(p);
        while (n % p == 0) n /= p;
    }
    if (n > 1) out.push_back(n);
    return out;
}

[[nodiscard]] inline std::size_t primitive_root_prime(std::size_t p) {
    if (p == 2) return 1;
    if (!is_prime(p)) throw std::invalid_argument("primitive root requires prime N");
    const auto factors = prime_factors(p - 1);
    for (std::size_t g = 2; g < p; ++g) {
        bool ok = true;
        for (const auto q : factors) {
            if (pow_mod(g, (p - 1) / q, p) == 1) { ok = false; break; }
        }
        if (ok) return g;
    }
    throw std::runtime_error("failed to find primitive root");
}

[[nodiscard]] inline std::size_t modular_inverse(std::size_t a, std::size_t m) {
    if (m == 0) throw std::invalid_argument("modular inverse modulus must be nonzero");
    // Extended Euclid using signed wide-enough arithmetic for ordinary size_t domains.
    using S = std::int64_t;
    if (a > static_cast<std::size_t>(std::numeric_limits<S>::max()) ||
        m > static_cast<std::size_t>(std::numeric_limits<S>::max()))
        throw std::overflow_error("modular inverse input exceeds supported signed range");
    S t = 0, new_t = 1;
    S r = static_cast<S>(m), new_r = static_cast<S>(a % m);
    while (new_r != 0) {
        const S q = r / new_r;
        std::tie(t, new_t) = std::pair<S, S>{new_t, t - q * new_t};
        std::tie(r, new_r) = std::pair<S, S>{new_r, r - q * new_r};
    }
    if (r != 1) throw std::invalid_argument("modular inverse does not exist");
    if (t < 0) t += static_cast<S>(m);
    return static_cast<std::size_t>(t);
}

[[nodiscard]] inline std::pair<std::size_t, std::size_t> coprime_factor_split(std::size_t n) noexcept {
    if (n < 6) return {0, 0};
    std::pair<std::size_t, std::size_t> best{0, 0};
    std::size_t best_gap = n;
    for (std::size_t a = 2; a <= n / a; ++a) {
        if (n % a != 0) continue;
        const auto b = n / a;
        if (std::gcd(a, b) != 1) continue;
        const auto gap = b - a;
        if (gap < best_gap) { best = {a, b}; best_gap = gap; }
    }
    return best;
}

[[nodiscard]] inline std::string_view name(Algo a) noexcept {
    switch (a) {
        case Algo::Auto: return "auto";
        case Algo::Dft: return "dft";
        case Algo::Radix2: return "radix2-iterative";
        case Algo::Recursive: return "radix2-recursive";
        case Algo::Stockham: return "stockham-radix2";
        case Algo::Radix4: return "radix4";
        case Algo::SplitRadix: return "split-radix";
        case Algo::ModifiedSplitRadix: return "modified-split-radix";
        case Algo::Mixed: return "mixed-radix";
        case Algo::GoodThomas: return "good-thomas";
        case Algo::Rader: return "rader";
        case Algo::Bluestein: return "bluestein";
    }
    return "unknown";
}

[[nodiscard]] inline Algo parse_algo(std::string_view s) {
    for (const auto a : all_algos) if (name(a) == s) return a;
    if (s == "radix2") return Algo::Radix2;
    if (s == "recursive") return Algo::Recursive;
    if (s == "stockham") return Algo::Stockham;
    if (s == "split") return Algo::SplitRadix;
    if (s == "modified-split") return Algo::ModifiedSplitRadix;
    if (s == "mixed") return Algo::Mixed;
    if (s == "pfa") return Algo::GoodThomas;
    if (s == "chirpz") return Algo::Bluestein;
    throw std::invalid_argument("unknown FFT algorithm");
}

[[nodiscard]] inline bool supports(Algo a, std::size_t n) noexcept {
    switch (a) {
        case Algo::Radix2:
        case Algo::Recursive:
        case Algo::Stockham:
        case Algo::Radix4:
        case Algo::SplitRadix:
        case Algo::ModifiedSplitRadix:
            return n == 0 || pow2(n);
        case Algo::GoodThomas: {
            const auto [x, y] = coprime_factor_split(n);
            return n <= 1 || (x > 1 && y > 1);
        }
        case Algo::Rader: return n <= 2 || is_prime(n);
        default: return true;
    }
}


} // namespace fftlab
