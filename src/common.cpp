#include "fftlab/fft.hpp"
#include "internal.hpp"

#include <algorithm>
#include <cmath>
#include <random>
#include <stdexcept>

namespace fftlab {
bool pow2(std::size_t n) { return n && !(n & (n - 1)); }
std::size_t ilog2(std::size_t n) {
    if (!pow2(n)) throw std::invalid_argument("ilog2 requires a power of two");
    std::size_t r = 0;
    while (n > 1) { n >>= 1; ++r; }
    return r;
}
std::size_t next_pow2(std::size_t n) {
    if (n <= 1) return 1;
    const auto high = std::size_t{1} << (std::numeric_limits<std::size_t>::digits - 1);
    if (n > high) throw std::overflow_error("next power of two overflows size_t");
    --n;
    for (std::size_t s = 1; s < std::numeric_limits<std::size_t>::digits; s <<= 1) n |= n >> s;
    return n + 1;
}
Complex root(double a) { return {std::cos(a), std::sin(a)}; }

bool is_prime(std::size_t n) {
    if (n < 2) return false;
    if (!(n % 2)) return n == 2;
    for (std::size_t d = 3; d <= n / d; d += 2) if (!(n % d)) return false;
    return true;
}
std::vector<std::size_t> prime_factors(std::size_t n) {
    std::vector<std::size_t> out;
    for (std::size_t p = 2; p <= n / p; ++p) {
        if (n % p) continue;
        out.push_back(p);
        while (!(n % p)) n /= p;
    }
    if (n > 1) out.push_back(n);
    return out;
}
std::size_t mul_mod(std::size_t a, std::size_t b, std::size_t m) {
    std::size_t r = 0;
    while (b) {
        if (b & 1) r = a >= m - r ? a - (m - r) : a + r;
        b >>= 1;
        if (b) a = a >= m - a ? a - (m - a) : a + a;
    }
    return r;
}
std::size_t pow_mod(std::size_t a, std::size_t e, std::size_t m) {
    std::size_t r = 1 % m;
    while (e) {
        if (e & 1) r = mul_mod(r, a, m);
        e >>= 1;
        if (e) a = mul_mod(a, a, m);
    }
    return r;
}
std::size_t primitive_root_prime(std::size_t p) {
    if (p == 2) return 1;
    if (!is_prime(p)) throw std::invalid_argument("primitive root requires prime N");
    const auto factors = prime_factors(p - 1);
    for (std::size_t g = 2; g < p; ++g) {
        bool ok = true;
        for (auto q : factors) if (pow_mod(g, (p - 1) / q, p) == 1) { ok = false; break; }
        if (ok) return g;
    }
    throw std::runtime_error("failed to find primitive root");
}

Vector dft(const Vector& x, bool inv) {
    const auto n = x.size(); Vector y(n);
    if (!n) return y;
    const double sign = inv ? 1.0 : -1.0, scale = inv ? 1.0 / double(n) : 1.0;
    for (std::size_t k = 0; k < n; ++k) {
        Complex s{};
        for (std::size_t t = 0; t < n; ++t)
            s += x[t] * root(sign * 2.0 * pi * double(k) * double(t) / double(n));
        y[k] = s * scale;
    }
    return y;
}

bool smooth235(std::size_t n) { if (!n) return false; for (auto p : {2u, 3u, 5u}) while (!(n % p)) n /= p; return n == 1; }
std::string_view name(Algo a) {
    switch (a) {
        case Algo::Auto: return "auto"; case Algo::Dft: return "dft"; case Algo::Radix2: return "radix2-iterative";
        case Algo::Recursive: return "radix2-recursive"; case Algo::Stockham: return "stockham-radix2"; case Algo::Radix4: return "radix4";
        case Algo::SplitRadix: return "split-radix"; case Algo::Mixed: return "mixed-radix"; case Algo::Rader: return "rader"; case Algo::Bluestein: return "bluestein";
    }
    return "?";
}
Algo parse_algo(std::string_view s) {
    if (s == "auto") return Algo::Auto;
    if (s == "dft") return Algo::Dft;
    if (s == "radix2" || s == "radix2-iterative") return Algo::Radix2;
    if (s == "recursive" || s == "radix2-recursive") return Algo::Recursive;
    if (s == "stockham" || s == "stockham-radix2") return Algo::Stockham;
    if (s == "radix4") return Algo::Radix4;
    if (s == "split" || s == "split-radix") return Algo::SplitRadix;
    if (s == "mixed" || s == "mixed-radix") return Algo::Mixed;
    if (s == "rader") return Algo::Rader;
    if (s == "bluestein" || s == "chirpz") return Algo::Bluestein;
    throw std::invalid_argument("unknown algorithm");
}
bool supports(Algo a, std::size_t n) {
    if (a == Algo::Radix2 || a == Algo::Recursive || a == Algo::Stockham || a == Algo::Radix4 || a == Algo::SplitRadix) return !n || pow2(n);
    if (a == Algo::Rader) return n <= 2 || is_prime(n);
    return true;
}
Vector transform(const Vector& x, Algo a, bool inv) {
    if (a == Algo::Auto) {
        if (x.size() <= 1) return x;
        if (pow2(x.size())) return radix2(x, inv);
        if (smooth235(x.size())) return mixed(x, inv);
        if (is_prime(x.size()) && x.size() >= 17) return rader(x, inv);
        return bluestein(x, inv);
    }
    if (a == Algo::Dft) return dft(x, inv);
    if (a == Algo::Radix2) return radix2(x, inv);
    if (a == Algo::Recursive) return recursive(x, inv);
    if (a == Algo::Stockham) return stockham(x, inv);
    if (a == Algo::Radix4) return radix4(x, inv);
    if (a == Algo::SplitRadix) return split_radix(x, inv);
    if (a == Algo::Mixed) return mixed(x, inv);
    if (a == Algo::Rader) return rader(x, inv);
    return bluestein(x, inv);
}

std::string_view signal_name(SignalKind s) {
    switch (s) { case SignalKind::Random: return "random"; case SignalKind::Tones: return "tones"; case SignalKind::Impulse: return "impulse"; case SignalKind::Alternating: return "alternating"; case SignalKind::DynamicRange: return "dynamic-range"; }
    return "?";
}
SignalKind parse_signal(std::string_view s) {
    for (auto k : signal_kinds) if (signal_name(k) == s) return k;
    throw std::invalid_argument("unknown signal kind");
}
Vector signal(std::size_t n, SignalKind kind, std::uint64_t seed) {
    Vector x(n); if (!n) return x;
    std::mt19937_64 g(seed ^ (0x9E3779B97F4A7C15ULL * n)); std::uniform_real_distribution<double> d(-0.5, 0.5);
    if (kind == SignalKind::Random) { for (auto& z : x) z = {d(g), d(g)}; return x; }
    if (kind == SignalKind::Impulse) { x[n / 3] = {1.0, -0.25}; return x; }
    if (kind == SignalKind::Alternating) { for (std::size_t i = 0; i < n; ++i) x[i] = i & 1 ? Complex{-1.0, 1e-12} : Complex{1.0, -1e-12}; return x; }
    if (kind == SignalKind::DynamicRange) {
        for (std::size_t i = 0; i < n; ++i) {
            const int decade = int(i % 25) - 12; const double scale = std::pow(10.0, double(decade));
            x[i] = scale * Complex{d(g), d(g)};
        }
        return x;
    }
    for (std::size_t i = 0; i < n; ++i) {
        const double t = double(i) / double(n);
        x[i] = {.9 * std::sin(2 * pi * 3 * t) + .35 * std::cos(2 * pi * 11 * t) + .05 * d(g), .4 * std::sin(2 * pi * 5 * t) + .05 * d(g)};
    }
    return x;
}

}
