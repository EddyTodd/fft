#include "fftlab/fft.hpp"
#include "internal.hpp"

#include <array>
#include <stdexcept>
#include <utility>

namespace fftlab {
void radix2_inplace(Vector& a, bool inv) {
    const auto n = a.size();
    if (!n) return;
    if (!pow2(n)) throw std::invalid_argument("radix2 requires power-of-two N");
    for (std::size_t i = 1, j = 0; i < n; ++i) {
        std::size_t bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) std::swap(a[i], a[j]);
    }
    const double sign = inv ? 1.0 : -1.0;
    for (std::size_t len = 2; len <= n;) {
        const Complex step = root(sign * 2.0 * pi / double(len));
        for (std::size_t i = 0; i < n; i += len) {
            Complex w{1, 0};
            for (std::size_t j = 0; j < len / 2; ++j) {
                const Complex u = a[i + j], v = a[i + j + len / 2] * w;
                a[i + j] = u + v;
                a[i + j + len / 2] = u - v;
                w *= step;
            }
        }
        if (len == n) break;
        len <<= 1;
    }
    if (inv) for (auto& z : a) z /= double(n);
}
Vector radix2(const Vector& x, bool inv) { Vector y = x; radix2_inplace(y, inv); return y; }

void recursive_core(Vector& a, bool inv) {
    const auto n = a.size(); if (n <= 1) return;
    Vector e(n / 2), o(n / 2);
    for (std::size_t i = 0; i < n / 2; ++i) { e[i] = a[2 * i]; o[i] = a[2 * i + 1]; }
    recursive_core(e, inv); recursive_core(o, inv);
    Complex w{1, 0}, step = root((inv ? 1.0 : -1.0) * 2.0 * pi / double(n));
    for (std::size_t k = 0; k < n / 2; ++k) {
        const auto t = w * o[k]; a[k] = e[k] + t; a[k + n / 2] = e[k] - t; w *= step;
    }
}
Vector recursive(const Vector& x, bool inv) {
    if (!x.empty() && !pow2(x.size())) throw std::invalid_argument("recursive radix2 requires power-of-two N");
    Vector y = x; recursive_core(y, inv); if (inv && !y.empty()) for (auto& z : y) z /= double(y.size()); return y;
}

Vector stockham(const Vector& x, bool inv) {
    const auto n = x.size();
    if (!n) return {};
    if (!pow2(n)) throw std::invalid_argument("Stockham requires power-of-two N");
    Vector src = x, dst(n);
    const double sign = inv ? 1.0 : -1.0;
    for (std::size_t m = 1; m < n; m <<= 1) {
        for (std::size_t j = 0; j < n / 2; ++j) {
            const auto k = j & (m - 1);
            const auto t = src[j + n / 2] * root(sign * pi * double(k) / double(m));
            const auto out = 2 * j - k;
            dst[out] = src[j] + t;
            dst[out + m] = src[j] - t;
        }
        src.swap(dst);
    }
    if (inv) for (auto& z : src) z /= double(n);
    return src;
}

void radix4_core(Vector& a, bool inv) {
    const auto n = a.size();
    if (n <= 1) return;
    if (n == 2) {
        const auto u = a[0], v = a[1]; a[0] = u + v; a[1] = u - v; return;
    }
    const auto m = n / 4;
    std::array<Vector, 4> sub{Vector(m), Vector(m), Vector(m), Vector(m)};
    for (std::size_t q = 0; q < 4; ++q) for (std::size_t j = 0; j < m; ++j) sub[q][j] = a[4 * j + q];
    for (auto& s : sub) radix4_core(s, inv);
    const double sign = inv ? 1.0 : -1.0;
    const Complex I{0.0, sign};
    for (std::size_t k = 0; k < m; ++k) {
        const auto b = sub[1][k] * root(sign * 2.0 * pi * double(k) / double(n));
        const auto c = sub[2][k] * root(sign * 4.0 * pi * double(k) / double(n));
        const auto d = sub[3][k] * root(sign * 6.0 * pi * double(k) / double(n));
        const auto t0 = sub[0][k] + c;
        const auto t1 = sub[0][k] - c;
        const auto t2 = b + d;
        const auto t3 = (b - d) * I;
        a[k] = t0 + t2;
        a[k + m] = t1 + t3;
        a[k + 2 * m] = t0 - t2;
        a[k + 3 * m] = t1 - t3;
    }
}
Vector radix4(const Vector& x, bool inv) {
    if (!x.empty() && !pow2(x.size())) throw std::invalid_argument("radix4 requires power-of-two N");
    Vector y = x; radix4_core(y, inv); if (inv && !y.empty()) for (auto& z : y) z /= double(y.size()); return y;
}

void split_core(Vector& a, bool inv) {
    const auto n = a.size();
    if (n <= 1) return;
    if (n == 2) { const auto u = a[0], v = a[1]; a[0] = u + v; a[1] = u - v; return; }
    Vector e(n / 2), o1(n / 4), o3(n / 4);
    for (std::size_t j = 0; j < n / 2; ++j) e[j] = a[2 * j];
    for (std::size_t j = 0; j < n / 4; ++j) { o1[j] = a[4 * j + 1]; o3[j] = a[4 * j + 3]; }
    split_core(e, inv); split_core(o1, inv); split_core(o3, inv);
    const double sign = inv ? 1.0 : -1.0;
    const Complex I{0.0, sign};
    for (std::size_t k = 0; k < n / 4; ++k) {
        const auto t1 = o1[k] * root(sign * 2.0 * pi * double(k) / double(n));
        const auto t2 = o3[k] * root(sign * 6.0 * pi * double(k) / double(n));
        const auto s = t1 + t2, d = t1 - t2;
        a[k] = e[k] + s;
        a[k + n / 2] = e[k] - s;
        a[k + n / 4] = e[k + n / 4] + I * d;
        a[k + 3 * n / 4] = e[k + n / 4] - I * d;
    }
}
Vector split_radix(const Vector& x, bool inv) {
    if (!x.empty() && !pow2(x.size())) throw std::invalid_argument("split-radix requires power-of-two N");
    Vector y = x; split_core(y, inv); if (inv && !y.empty()) for (auto& z : y) z /= double(y.size()); return y;
}

}
