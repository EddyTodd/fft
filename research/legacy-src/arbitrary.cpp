#include "fftlab/fft.hpp"
#include "internal.hpp"

#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <utility>

namespace fftlab {
std::size_t factor(std::size_t n) {
    if (!(n % 2)) return 2;
    for (std::size_t p = 3; p <= n / p; p += 2) if (!(n % p)) return p;
    return n;
}
Vector mixed_core(const Vector& x, bool inv) {
    const auto n = x.size(); if (n <= 1) return x;
    const auto r = factor(n);
    if (r == n) {
        Vector y(n); const double sign = inv ? 1.0 : -1.0;
        for (std::size_t k = 0; k < n; ++k) for (std::size_t t = 0; t < n; ++t)
            y[k] += x[t] * root(sign * 2.0 * pi * double(k) * double(t) / double(n));
        return y;
    }
    const auto m = n / r; std::vector<Vector> sub(r, Vector(m));
    for (std::size_t q = 0; q < r; ++q) {
        for (std::size_t j = 0; j < m; ++j) sub[q][j] = x[r * j + q];
        sub[q] = mixed_core(sub[q], inv);
    }
    Vector y(n); const double sign = inv ? 1.0 : -1.0;
    for (std::size_t k0 = 0; k0 < m; ++k0) for (std::size_t k1 = 0; k1 < r; ++k1) {
        const auto k = k0 + m * k1; Complex w{1, 0}, step = root(sign * 2.0 * pi * double(k) / double(n));
        for (std::size_t q = 0; q < r; ++q) { y[k] += sub[q][k0] * w; w *= step; }
    }
    return y;
}
Vector mixed(const Vector& x, bool inv) {
    Vector y = mixed_core(x, inv); if (inv && !y.empty()) for (auto& z : y) z /= double(y.size()); return y;
}

Vector bluestein(const Vector& x, bool inv) {
    const auto n = x.size(); if (n <= 1) return x;
    if (n > (std::numeric_limits<std::size_t>::max() / 2) + 1) throw std::length_error("Bluestein workspace overflow");
    const auto m = next_pow2(2 * n - 1); Vector a(m), b(m); const double sign = inv ? 1.0 : -1.0;
    for (std::size_t k = 0; k < n; ++k) {
        const long double kd = k, period = 2.0L * n, phase = std::fmod(kd * kd, period) / n;
        const Complex c = root(sign * pi * double(phase)); a[k] = x[k] * c; b[k] = std::conj(c); if (k) b[m - k] = std::conj(c);
    }
    radix2_inplace(a); radix2_inplace(b); for (std::size_t i = 0; i < m; ++i) a[i] *= b[i]; radix2_inplace(a, true);
    Vector y(n);
    for (std::size_t k = 0; k < n; ++k) {
        const long double kd = k, period = 2.0L * n, phase = std::fmod(kd * kd, period) / n;
        y[k] = a[k] * root(sign * pi * double(phase)); if (inv) y[k] /= double(n);
    }
    return y;
}

Vector circular_convolution(Vector a, Vector b) {
    if (a.size() != b.size()) throw std::invalid_argument("circular convolution size mismatch");
    const auto n = a.size(); if (!n) return {};
    const auto m = next_pow2(2 * n - 1); a.resize(m); b.resize(m);
    radix2_inplace(a); radix2_inplace(b);
    for (std::size_t i = 0; i < m; ++i) a[i] *= b[i];
    radix2_inplace(a, true);
    Vector c(n);
    for (std::size_t i = 0; i < 2 * n - 1; ++i) c[i % n] += a[i];
    return c;
}
Vector rader(const Vector& x, bool inv) {
    const auto n = x.size();
    if (n <= 2) return dft(x, inv);
    if (!is_prime(n)) throw std::invalid_argument("Rader requires prime N");
    const auto g = primitive_root_prime(n), l = n - 1;
    Vector ar(l), b(l); std::vector<std::size_t> gp(l);
    gp[0] = 1;
    for (std::size_t q = 1; q < l; ++q) gp[q] = mul_mod(gp[q - 1], g, n);
    const double sign = inv ? 1.0 : -1.0;
    for (std::size_t q = 0; q < l; ++q) {
        const auto negq = (l - q) % l;
        ar[q] = x[gp[negq]];
        b[q] = root(sign * 2.0 * pi * double(gp[q]) / double(n));
    }
    const auto c = circular_convolution(std::move(ar), std::move(b));
    Vector y(n); y[0] = std::accumulate(x.begin(), x.end(), Complex{});
    for (std::size_t m = 0; m < l; ++m) y[gp[m]] = x[0] + c[m];
    if (inv) for (auto& z : y) z /= double(n);
    return y;
}

}
