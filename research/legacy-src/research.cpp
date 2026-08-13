#include "fftlab/fft.hpp"
#include "internal.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>

namespace fftlab {
using LComplex = std::complex<long double>;
using LVector = std::vector<LComplex>;
LVector oracle(const Vector& x, bool inv = false) {
    const auto n = x.size(); LVector y(n); if (!n) return y;
    const long double sign = inv ? 1.0L : -1.0L, scale = inv ? 1.0L / n : 1.0L;
    for (std::size_t k = 0; k < n; ++k) for (std::size_t t = 0; t < n; ++t) {
        const long double a = sign * 2 * lpi * k * t / n;
        y[k] += LComplex{x[t].real(), x[t].imag()} * LComplex{std::cos(a), std::sin(a)};
    }
    if (inv) for (auto& z : y) z *= scale;
    return y;
}
LVector widen(const Vector& x) { LVector y; y.reserve(x.size()); for (auto z : x) y.emplace_back(z.real(), z.imag()); return y; }
NormErr norm_error(const LVector& got, const LVector& ref) {
    if (got.size() != ref.size()) throw std::invalid_argument("error vector size mismatch");
    long double e1 = 0, e2 = 0, ei = 0, r1 = 0, r2 = 0, ri = 0;
    for (std::size_t i = 0; i < got.size(); ++i) {
        const auto e = std::abs(got[i] - ref[i]), r = std::abs(ref[i]);
        e1 += e; e2 += e * e; ei = std::max(ei, e); r1 += r; r2 += r * r; ri = std::max(ri, r);
    }
    return {double(r1 != 0 ? e1 / r1 : e1), double(r2 != 0 ? std::sqrt(e2 / r2) : std::sqrt(e2)), double(ri != 0 ? ei / ri : ei)};
}
Accuracy accuracy(const Vector& x, Algo a, SignalKind sig) {
    const auto exact = oracle(x);
    const auto y = transform(x, a);
    const auto back_exact = oracle(y, true);
    double rt = 0; const auto rtvec = transform(y, a, true);
    for (std::size_t i = 0; i < x.size(); ++i) rt = std::max(rt, std::abs(rtvec[i] - x[i]));
    return {a, x.size(), sig, norm_error(widen(y), exact), norm_error(back_exact, widen(x)), rt};
}

Model dft_model(std::size_t n) { const auto nn = static_cast<long double>(n); return {nn * static_cast<long double>(n ? n - 1 : 0), nn * nn, static_cast<long double>(n), "direct DFT structural count"}; }
Model radix2_model(std::size_t n) {
    if (!n) return {};
    const auto l = ilog2(n); return {static_cast<long double>(n) * static_cast<long double>(l), static_cast<long double>(n) * static_cast<long double>(l) / 2, 0, "canonical radix-2 butterflies; excludes twiddle-generation arithmetic"};
}
Model radix4_model(std::size_t n) {
    if (n <= 1) return {};
    if (n == 2) return {2, 0, 0, "radix-2 base case"};
    auto s = radix4_model(n / 4); return {4 * s.adds + 2 * static_cast<long double>(n), 4 * s.muls + 0.75L * static_cast<long double>(n), static_cast<long double>(n), "optimized radix-4 butterfly; +/-i treated as trivial"};
}
Model split_model(std::size_t n) {
    if (n <= 1) return {};
    if (n == 2) return {2, 0, 0, "radix-2 base case"};
    auto h = split_model(n / 2), q = split_model(n / 4);
    return {h.adds + 2 * q.adds + 1.5L * static_cast<long double>(n), h.muls + 2 * q.muls + 0.5L * static_cast<long double>(n), static_cast<long double>(n), "split-radix recurrence; +/-i treated as trivial"};
}
Model mixed_model(std::size_t n) {
    if (n <= 1) return {};
    const auto r = factor(n); if (r == n) return dft_model(n);
    auto s = mixed_model(n / r);
    return {static_cast<long double>(r) * s.adds + static_cast<long double>(n) * static_cast<long double>(r - 1), static_cast<long double>(r) * s.muls + static_cast<long double>(n) * static_cast<long double>(r), static_cast<long double>(2 * n), "structural count for this pedagogical mixed-radix implementation"};
}
Model model(Algo a, std::size_t n) {
    if (!supports(a, n)) throw std::invalid_argument("algorithm does not support N");
    if (a == Algo::Auto) a = pow2(n) ? Algo::Radix2 : (smooth235(n) ? Algo::Mixed : (is_prime(n) && n >= 17 ? Algo::Rader : Algo::Bluestein));
    if (a == Algo::Dft) return dft_model(n);
    if (a == Algo::Radix2 || a == Algo::Recursive || a == Algo::Stockham) { auto m = radix2_model(n); if (a == Algo::Stockham) m.workspace = n; if (a == Algo::Recursive) m.workspace = 2 * n; return m; }
    if (a == Algo::Radix4) return radix4_model(n);
    if (a == Algo::SplitRadix) return split_model(n);
    if (a == Algo::Mixed) return mixed_model(n);
    if (a == Algo::Bluestein) {
        if (n <= 1) return {};
        const auto m = next_pow2(2 * n - 1); auto r = radix2_model(m);
        return {3 * r.adds, 3 * r.muls + static_cast<long double>(m) + 2 * static_cast<long double>(n), static_cast<long double>(2 * m + n), "three radix-2 FFTs + pointwise convolution + chirps"};
    }
    if (n <= 2) return dft_model(n);
    const auto l = n - 1, m = next_pow2(2 * l - 1); auto r = radix2_model(m);
    return {3 * r.adds + static_cast<long double>(l), 3 * r.muls + static_cast<long double>(m), static_cast<long double>(2 * m + 2 * l + n), "Rader reduction to cyclic convolution; permutation/root setup excluded"};
}

volatile double sink = 0;
double pct(const std::vector<double>& v, double p) {
    const double pos = p * double(v.size() - 1); const auto lo = std::size_t(std::floor(pos)), hi = std::size_t(std::ceil(pos));
    return lo == hi ? v[lo] : v[lo] * (double(hi) - pos) + v[hi] * (pos - double(lo));
}
std::pair<double, double> bootstrap_median_ci(const std::vector<double>& values, std::uint64_t seed, std::size_t reps = 2000) {
    std::mt19937_64 rng(seed); std::uniform_int_distribution<std::size_t> pick(0, values.size() - 1); std::vector<double> boots; boots.reserve(reps); std::vector<double> sample(values.size());
    for (std::size_t b = 0; b < reps; ++b) { for (auto& x : sample) x = values[pick(rng)]; std::sort(sample.begin(), sample.end()); boots.push_back(pct(sample, .5)); }
    std::sort(boots.begin(), boots.end()); return {pct(boots, .025), pct(boots, .975)};
}
std::size_t calibrate(const Vector& x, Algo a, double ms) {
    using C = std::chrono::steady_clock; std::size_t it = 1;
    while (it < (1u << 20)) {
        const auto s = C::now(); double c = 0;
        for (std::size_t i = 0; i < it; ++i) { auto y = transform(x, a); if (!y.empty()) c += y[i % y.size()].real(); }
        sink = c; const double e = std::chrono::duration<double, std::milli>(C::now() - s).count(); if (e >= ms * .5) break;
        it = std::min<std::size_t>(1u << 20, std::max(it + 1, std::size_t(double(it) * std::clamp(ms / std::max(e, 1e-9), 2.0, 16.0))));
    }
    return it;
}
Stats bench(const Vector& x, Algo a, std::size_t samples, std::size_t warm, double target) {
    if (!supports(a, x.size())) throw std::invalid_argument("unsupported size");
    if (samples < 5 || target <= 0) throw std::invalid_argument("bad benchmark parameters");
    for (std::size_t i = 0; i < warm; ++i) { auto y = transform(x, a); if (!y.empty()) sink = y[0].real(); }
    const auto it = calibrate(x, a, target); std::vector<double> raw; raw.reserve(samples); using C = std::chrono::steady_clock;
    for (std::size_t s = 0; s < samples; ++s) {
        const auto st = C::now(); double c = 0;
        for (std::size_t i = 0; i < it; ++i) { auto y = transform(x, a); if (!y.empty()) c += y[(i + s) % y.size()].real(); }
        sink = c; raw.push_back(std::chrono::duration<double, std::nano>(C::now() - st).count() / double(it));
    }
    auto sorted = raw; std::sort(sorted.begin(), sorted.end()); const double mean = std::accumulate(sorted.begin(), sorted.end(), 0.0) / static_cast<double>(sorted.size()); double var = 0;
    for (double z : sorted) var += (z - mean) * (z - mean);
    var /= static_cast<double>(sorted.size() - 1);
    const double med = pct(sorted, .5);
    std::vector<double> dev; dev.reserve(sorted.size()); for (double z : sorted) dev.push_back(std::abs(z - med)); std::sort(dev.begin(), dev.end());
    const auto ci = bootstrap_median_ci(sorted, 0xB00757A9ULL ^ (std::uint64_t(x.size()) << 8) ^ std::uint64_t(a));
    const double scaled = x.size() > 1 ? 5.0 * double(x.size()) * std::log2(double(x.size())) / med : 0.0;
    return {x.size(), samples, it, a, sorted.front(), pct(sorted, .05), med, mean, pct(sorted, .95), sorted.back(), std::sqrt(var), pct(dev, .5), ci.first, ci.second, scaled, std::move(raw)};
}

std::vector<Algo> suite(std::size_t n) {
    std::vector<Algo> v{Algo::Auto, Algo::Mixed, Algo::Bluestein}; if (is_prime(n)) v.push_back(Algo::Rader);
    if (pow2(n)) { v.push_back(Algo::Radix2); v.push_back(Algo::Recursive); v.push_back(Algo::Stockham); v.push_back(Algo::Radix4); v.push_back(Algo::SplitRadix); }
    if (n <= 2048) v.push_back(Algo::Dft);
    return v;
}

void tests() {
    std::size_t checks = 0; auto req = [&](bool c, std::string_view what) { ++checks; if (!c) throw std::runtime_error(std::string("self-test failure: ") + std::string(what)); };
    req(next_pow2(3) == 4, "next_pow2"); req(is_prime(1009), "prime"); req(!is_prime(1000), "composite"); req(primitive_root_prime(17) > 1, "primitive-root");
    for (std::size_t n : {1u, 2u, 3u, 5u, 6u, 7u, 8u, 10u, 12u, 16u, 17u, 25u, 31u, 32u, 64u}) {
        for (auto sk : signal_kinds) {
            const auto x = signal(n, sk), ref = dft(x);
            for (auto a : all_algos) {
                if (a == Algo::Dft || !supports(a, n)) continue;
                const auto y = transform(x, a); req(y.size() == n, "size");
                for (std::size_t i = 0; i < n; ++i) req(std::abs(y[i] - ref[i]) < 1e-8 * (1 + std::abs(ref[i])), name(a));
                const auto z = transform(y, a, true);
                double xscale = 0.0; for (auto value : x) xscale = std::max(xscale, std::abs(value));
                for (std::size_t i = 0; i < n; ++i) req(std::abs(z[i] - x[i]) < 2e-8 * (1 + xscale), std::string("inverse ") + std::string(name(a)) + " N=" + std::to_string(n));
            }
        }
    }
    bool threw = false; try { (void)rader(signal(15)); } catch (const std::invalid_argument&) { threw = true; } req(threw, "Rader rejects composite");
    threw = false; try { (void)stockham(signal(6)); } catch (const std::invalid_argument&) { threw = true; } req(threw, "Stockham rejects composite");
    std::cout << "PASS: " << checks << " checks\n";
}
}
