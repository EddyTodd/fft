#pragma once

#include <array>
#include <complex>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string_view>
#include <tuple>
#include <utility>
#include <vector>

namespace fftlab {

template <class T>
concept FftScalar = std::same_as<T, float> || std::same_as<T, double>;

template <FftScalar T>
using ComplexT = std::complex<T>;

template <FftScalar T>
using VectorT = std::vector<ComplexT<T>>;

using Complex32 = ComplexT<float>;
using Complex64 = ComplexT<double>;
using Vector32 = VectorT<float>;
using Vector64 = VectorT<double>;

enum class Direction {
    Forward,
    Inverse,
};

enum class Algorithm {
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
    Bluestein,
};

inline constexpr std::array<Algorithm, 12> all_algorithms{
    Algorithm::Auto,
    Algorithm::Dft,
    Algorithm::Radix2,
    Algorithm::Recursive,
    Algorithm::Stockham,
    Algorithm::Radix4,
    Algorithm::SplitRadix,
    Algorithm::ModifiedSplitRadix,
    Algorithm::Mixed,
    Algorithm::GoodThomas,
    Algorithm::Rader,
    Algorithm::Bluestein,
};

[[nodiscard]] inline constexpr bool pow2(std::size_t n) noexcept {
    return n != 0 && (n & (n - 1)) == 0;
}

[[nodiscard]] inline bool is_prime(std::size_t n) noexcept {
    if (n < 2) {
        return false;
    }
    if ((n & 1U) == 0U) {
        return n == 2;
    }
    for (std::size_t divisor = 3; divisor <= n / divisor; divisor += 2) {
        if (n % divisor == 0) {
            return false;
        }
    }
    return true;
}

[[nodiscard]] inline std::size_t next_pow2(std::size_t n) {
    if (n <= 1) {
        return 1;
    }

    const auto high = std::size_t{1} << (std::numeric_limits<std::size_t>::digits - 1);
    if (n > high) {
        throw std::overflow_error("next power of two overflows size_t");
    }

    --n;
    for (std::size_t shift = 1; shift < std::numeric_limits<std::size_t>::digits; shift <<= 1) {
        n |= n >> shift;
    }
    return n + 1;
}

[[nodiscard]] inline std::size_t ilog2(std::size_t n) {
    if (!pow2(n)) {
        throw std::invalid_argument("ilog2 requires a power of two");
    }

    std::size_t result = 0;
    while (n > 1) {
        n >>= 1;
        ++result;
    }
    return result;
}

[[nodiscard]] inline std::size_t mul_mod(std::size_t a, std::size_t b,
                                         std::size_t modulus) noexcept {
    std::size_t result = 0;
    while (b != 0) {
        if ((b & 1U) != 0U) {
            result = a >= modulus - result ? a - (modulus - result) : a + result;
        }
        b >>= 1;
        if (b != 0) {
            a = a >= modulus - a ? a - (modulus - a) : a + a;
        }
    }
    return result;
}

[[nodiscard]] inline std::size_t pow_mod(std::size_t base, std::size_t exponent,
                                         std::size_t modulus) noexcept {
    std::size_t result = 1 % modulus;
    while (exponent != 0) {
        if ((exponent & 1U) != 0U) {
            result = mul_mod(result, base, modulus);
        }
        exponent >>= 1;
        if (exponent != 0) {
            base = mul_mod(base, base, modulus);
        }
    }
    return result;
}

[[nodiscard]] inline std::vector<std::size_t> prime_factors(std::size_t n) {
    std::vector<std::size_t> factors;
    for (std::size_t factor = 2; factor <= n / factor; ++factor) {
        if (n % factor != 0) {
            continue;
        }
        factors.push_back(factor);
        while (n % factor == 0) {
            n /= factor;
        }
    }
    if (n > 1) {
        factors.push_back(n);
    }
    return factors;
}

[[nodiscard]] inline std::size_t primitive_root_prime(std::size_t prime) {
    if (prime == 2) {
        return 1;
    }
    if (!is_prime(prime)) {
        throw std::invalid_argument("primitive root requires prime N");
    }

    const auto factors = prime_factors(prime - 1);
    for (std::size_t candidate = 2; candidate < prime; ++candidate) {
        bool valid = true;
        for (const auto factor : factors) {
            if (pow_mod(candidate, (prime - 1) / factor, prime) == 1) {
                valid = false;
                break;
            }
        }
        if (valid) {
            return candidate;
        }
    }
    throw std::runtime_error("failed to find primitive root");
}

[[nodiscard]] inline std::size_t modular_inverse(std::size_t a, std::size_t modulus) {
    if (modulus == 0) {
        throw std::invalid_argument("modular inverse modulus must be nonzero");
    }

    using Signed = std::int64_t;
    if (a > static_cast<std::size_t>(std::numeric_limits<Signed>::max()) ||
        modulus > static_cast<std::size_t>(std::numeric_limits<Signed>::max())) {
        throw std::overflow_error("modular inverse input exceeds supported signed range");
    }

    Signed t = 0;
    Signed next_t = 1;
    Signed remainder = static_cast<Signed>(modulus);
    Signed next_remainder = static_cast<Signed>(a % modulus);
    while (next_remainder != 0) {
        const Signed quotient = remainder / next_remainder;
        std::tie(t, next_t) = std::pair<Signed, Signed>{next_t, t - quotient * next_t};
        std::tie(remainder, next_remainder) =
            std::pair<Signed, Signed>{next_remainder, remainder - quotient * next_remainder};
    }

    if (remainder != 1) {
        throw std::invalid_argument("modular inverse does not exist");
    }
    if (t < 0) {
        t += static_cast<Signed>(modulus);
    }
    return static_cast<std::size_t>(t);
}

[[nodiscard]] inline std::pair<std::size_t, std::size_t> coprime_factor_split(
    std::size_t n) noexcept {
    if (n < 6) {
        return {0, 0};
    }

    std::pair<std::size_t, std::size_t> best{0, 0};
    std::size_t best_gap = n;
    for (std::size_t first = 2; first <= n / first; ++first) {
        if (n % first != 0) {
            continue;
        }
        const auto second = n / first;
        if (std::gcd(first, second) != 1) {
            continue;
        }
        const auto gap = second - first;
        if (gap < best_gap) {
            best = {first, second};
            best_gap = gap;
        }
    }
    return best;
}

[[nodiscard]] inline std::string_view algorithm_name(Algorithm algorithm) noexcept {
    switch (algorithm) {
        case Algorithm::Auto:
            return "auto";
        case Algorithm::Dft:
            return "dft";
        case Algorithm::Radix2:
            return "radix2-iterative";
        case Algorithm::Recursive:
            return "radix2-recursive";
        case Algorithm::Stockham:
            return "stockham-radix2";
        case Algorithm::Radix4:
            return "radix4";
        case Algorithm::SplitRadix:
            return "split-radix";
        case Algorithm::ModifiedSplitRadix:
            return "modified-split-radix";
        case Algorithm::Mixed:
            return "mixed-radix";
        case Algorithm::GoodThomas:
            return "good-thomas";
        case Algorithm::Rader:
            return "rader";
        case Algorithm::Bluestein:
            return "bluestein";
    }
    return "unknown";
}

[[nodiscard]] inline Algorithm parse_algorithm(std::string_view name) {
    for (const auto algorithm : all_algorithms) {
        if (algorithm_name(algorithm) == name) {
            return algorithm;
        }
    }

    if (name == "radix2") {
        return Algorithm::Radix2;
    }
    if (name == "recursive") {
        return Algorithm::Recursive;
    }
    if (name == "stockham") {
        return Algorithm::Stockham;
    }
    if (name == "split") {
        return Algorithm::SplitRadix;
    }
    if (name == "modified-split") {
        return Algorithm::ModifiedSplitRadix;
    }
    if (name == "mixed") {
        return Algorithm::Mixed;
    }
    if (name == "pfa") {
        return Algorithm::GoodThomas;
    }
    if (name == "chirpz") {
        return Algorithm::Bluestein;
    }
    throw std::invalid_argument("unknown FFT algorithm");
}

[[nodiscard]] inline bool supports_algorithm(Algorithm algorithm, std::size_t n) noexcept {
    switch (algorithm) {
        case Algorithm::Radix2:
        case Algorithm::Recursive:
        case Algorithm::Stockham:
        case Algorithm::Radix4:
        case Algorithm::SplitRadix:
        case Algorithm::ModifiedSplitRadix:
            return n == 0 || pow2(n);
        case Algorithm::GoodThomas: {
            const auto [first, second] = coprime_factor_split(n);
            return n <= 1 || (first > 1 && second > 1);
        }
        case Algorithm::Rader:
            return n <= 2 || is_prime(n);
        default:
            return true;
    }
}

}  // namespace fftlab
