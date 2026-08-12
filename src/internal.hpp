#pragma once

#include "fftlab/fft.hpp"

#include <numbers>
#include <vector>

namespace fftlab {
inline constexpr double pi = std::numbers::pi_v<double>;
inline constexpr long double lpi = std::numbers::pi_v<long double>;
std::size_t ilog2(std::size_t n);
std::size_t next_pow2(std::size_t n);
Complex root(double a);
std::vector<std::size_t> prime_factors(std::size_t n);
std::size_t mul_mod(std::size_t a, std::size_t b, std::size_t m);
std::size_t pow_mod(std::size_t a, std::size_t e, std::size_t m);
std::size_t primitive_root_prime(std::size_t p);
bool smooth235(std::size_t n);
std::size_t factor(std::size_t n);
}
