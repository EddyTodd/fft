#include "fftlab/plan.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <numbers>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>

namespace fftlab {
namespace {
constexpr double pi = std::numbers::pi_v<double>;
Complex root(double angle) { return {std::cos(angle), std::sin(angle)}; }

double percentile(const std::vector<double>& sorted, double p) {
    const double pos = p * double(sorted.size() - 1);
    const auto lo = std::size_t(std::floor(pos)), hi = std::size_t(std::ceil(pos));
    return lo == hi ? sorted[lo] : sorted[lo] * (double(hi) - pos) + sorted[hi] * (pos - double(lo));
}
std::pair<double, double> bootstrap_median_ci(const std::vector<double>& values, std::uint64_t seed, std::size_t reps = 2000) {
    std::mt19937_64 rng(seed); std::uniform_int_distribution<std::size_t> pick(0, values.size() - 1);
    std::vector<double> boots; boots.reserve(reps); std::vector<double> sample(values.size());
    for (std::size_t b = 0; b < reps; ++b) {
        for (auto& x : sample) x = values[pick(rng)];
        std::sort(sample.begin(), sample.end()); boots.push_back(percentile(sample, .5));
    }
    std::sort(boots.begin(), boots.end()); return {percentile(boots, .025), percentile(boots, .975)};
}
PlanDistribution summarize(std::vector<double> raw, std::uint64_t seed) {
    auto sorted = raw; std::sort(sorted.begin(), sorted.end()); const double med = percentile(sorted, .5);
    std::vector<double> dev; dev.reserve(sorted.size()); for (double x : sorted) dev.push_back(std::abs(x - med)); std::sort(dev.begin(), dev.end());
    const auto ci = bootstrap_median_ci(sorted, seed);
    return {sorted.front(), percentile(sorted, .05), med, percentile(sorted, .95), sorted.back(), percentile(dev, .5), ci.first, ci.second, std::move(raw)};
}
volatile double sink = 0;
}

Radix2Plan::Radix2Plan(std::size_t n) : n_(n), bit_reverse_(n), twiddles_(n / 2) {
    if (!pow2(n)) throw std::invalid_argument("Radix2Plan requires power-of-two N");
    std::size_t bits = 0; for (auto q = n; q > 1; q >>= 1) ++bits;
    for (std::size_t i = 0; i < n; ++i) {
        auto x = i; std::size_t reversed = 0;
        for (std::size_t b = 0; b < bits; ++b) { reversed = (reversed << 1) | (x & 1); x >>= 1; }
        bit_reverse_[i] = reversed;
    }
    for (std::size_t k = 0; k < n / 2; ++k) twiddles_[k] = root(-2.0 * pi * double(k) / double(n));
}
void Radix2Plan::forward_inplace(Vector& data) const { execute(data, false); }
void Radix2Plan::inverse_inplace(Vector& data) const { execute(data, true); }
void Radix2Plan::execute(Vector& data, bool inverse) const {
    if (data.size() != n_) throw std::invalid_argument("Radix2Plan input size mismatch");
    for (std::size_t i = 0; i < n_; ++i) if (i < bit_reverse_[i]) std::swap(data[i], data[bit_reverse_[i]]);
    for (std::size_t len = 2; len <= n_;) {
        const auto stride = n_ / len;
        for (std::size_t base = 0; base < n_; base += len) {
            for (std::size_t j = 0; j < len / 2; ++j) {
                auto w = twiddles_[j * stride]; if (inverse) w = std::conj(w);
                const auto u = data[base + j], v = data[base + j + len / 2] * w;
                data[base + j] = u + v; data[base + j + len / 2] = u - v;
            }
        }
        if (len == n_) break;
        len <<= 1;
    }
    if (inverse) for (auto& z : data) z /= double(n_);
}

RealRadix2Plan::RealRadix2Plan(std::size_t n)
    : n_(n), half_(n / 2), half_plan_(n / 2), post_twiddles_(n / 2 + 1) {
    if (n < 2 || !pow2(n)) throw std::invalid_argument("RealRadix2Plan requires power-of-two N >= 2");
    for (std::size_t k = 0; k <= half_; ++k) post_twiddles_[k] = root(-2.0 * pi * double(k) / double(n_));
}
void RealRadix2Plan::forward(const RealVector& input, HalfSpectrum& output, Vector& scratch) const {
    if (input.size() != n_ || output.size() != half_ + 1 || scratch.size() != half_) throw std::invalid_argument("RealRadix2Plan buffer size mismatch");
    for (std::size_t j = 0; j < half_; ++j) scratch[j] = {input[2 * j], input[2 * j + 1]};
    half_plan_.forward_inplace(scratch);
    for (std::size_t k = 0; k <= half_; ++k) {
        const auto a = scratch[k % half_], b = std::conj(scratch[(half_ - k) % half_]);
        const auto even = 0.5 * (a + b), odd = Complex{0.0, -0.5} * (a - b);
        output[k] = even + post_twiddles_[k] * odd;
    }
}
void RealRadix2Plan::inverse(const HalfSpectrum& input, RealVector& output, Vector& scratch) const {
    if (input.size() != half_ + 1 || output.size() != n_ || scratch.size() != half_) throw std::invalid_argument("RealRadix2Plan buffer size mismatch");
    for (std::size_t k = 0; k < half_; ++k) {
        const auto xk = input[k], mirror = std::conj(input[k ? half_ - k : half_]);
        const auto even = 0.5 * (xk + mirror), odd = 0.5 * (xk - mirror) / post_twiddles_[k];
        scratch[k] = even + Complex{0.0, 1.0} * odd;
    }
    half_plan_.inverse_inplace(scratch);
    for (std::size_t j = 0; j < half_; ++j) { output[2 * j] = scratch[j].real(); output[2 * j + 1] = scratch[j].imag(); }
}

PlanBenchmark benchmark_plans(std::size_t n, std::size_t samples, std::size_t warmups, double target_ms, std::uint64_t seed) {
    if (n < 2 || !pow2(n)) throw std::invalid_argument("planned benchmark requires power-of-two N >= 2");
    if (samples < 5 || target_ms <= 0) throw std::invalid_argument("bad planned benchmark parameters");
    std::mt19937_64 rng(seed ^ n); std::uniform_real_distribution<double> dist(-0.5, 0.5);
    Vector complex_input(n); for (auto& z : complex_input) z = {dist(rng), dist(rng)};
    RealVector real_input(n); for (auto& x : real_input) x = dist(rng);
    Radix2Plan plan(n); RealRadix2Plan real_plan(n); Vector cbuf = complex_input, scratch(n / 2); HalfSpectrum spectrum(n / 2 + 1); RealVector rbuf = real_input;
    for (std::size_t i = 0; i < warmups; ++i) { plan.forward_inplace(cbuf); plan.inverse_inplace(cbuf); real_plan.forward(rbuf, spectrum, scratch); real_plan.inverse(spectrum, rbuf, scratch); }
    using Clock = std::chrono::steady_clock;
    std::size_t iters = 1;
    while (iters < (1u << 20)) {
        const auto start = Clock::now();
        for (std::size_t i = 0; i < iters; ++i) { plan.forward_inplace(cbuf); plan.inverse_inplace(cbuf); }
        const auto ms = std::chrono::duration<double, std::milli>(Clock::now() - start).count();
        if (ms >= target_ms) break;
        iters = std::min<std::size_t>(1u << 20, std::max(iters + 1, std::size_t(double(iters) * std::clamp(target_ms / std::max(ms, 1e-9), 2.0, 16.0))));
    }
    std::vector<double> complex_setup_raw, real_setup_raw, legacy_raw, planned_raw, real_raw;
    complex_setup_raw.reserve(samples); real_setup_raw.reserve(samples); legacy_raw.reserve(samples); planned_raw.reserve(samples); real_raw.reserve(samples);
    std::array<int, 5> modes{0, 1, 2, 3, 4};
    for (std::size_t s = 0; s < samples; ++s) {
        std::shuffle(modes.begin(), modes.end(), rng);
        for (int mode : modes) {
            if (mode == 0) {
                const auto start = Clock::now(); Radix2Plan fresh(n); const auto end = Clock::now(); sink = double(fresh.stored_twiddles());
                complex_setup_raw.push_back(std::chrono::duration<double, std::nano>(end - start).count());
            } else if (mode == 1) {
                const auto start = Clock::now(); RealRadix2Plan fresh(n); const auto end = Clock::now(); sink = double(fresh.spectrum_size());
                real_setup_raw.push_back(std::chrono::duration<double, std::nano>(end - start).count());
            } else if (mode == 2) {
                cbuf = complex_input; const auto start = Clock::now(); for (std::size_t i = 0; i < iters; ++i) { radix2_inplace(cbuf); radix2_inplace(cbuf, true); } const auto end = Clock::now(); sink = cbuf[0].real();
                legacy_raw.push_back(std::chrono::duration<double, std::nano>(end - start).count() / double(2 * iters));
            } else if (mode == 3) {
                cbuf = complex_input; const auto start = Clock::now(); for (std::size_t i = 0; i < iters; ++i) { plan.forward_inplace(cbuf); plan.inverse_inplace(cbuf); } const auto end = Clock::now(); sink = cbuf[0].real();
                planned_raw.push_back(std::chrono::duration<double, std::nano>(end - start).count() / double(2 * iters));
            } else {
                rbuf = real_input; const auto start = Clock::now(); for (std::size_t i = 0; i < iters; ++i) { real_plan.forward(rbuf, spectrum, scratch); real_plan.inverse(spectrum, rbuf, scratch); } const auto end = Clock::now(); sink = rbuf[0];
                real_raw.push_back(std::chrono::duration<double, std::nano>(end - start).count() / double(2 * iters));
            }
        }
    }
    auto complex_setup = summarize(std::move(complex_setup_raw), seed ^ 0x51E7ULL ^ n);
    auto real_setup = summarize(std::move(real_setup_raw), seed ^ 0x7E7ULL ^ n);
    auto legacy = summarize(std::move(legacy_raw), seed ^ 0x1E6AULL ^ n);
    auto planned = summarize(std::move(planned_raw), seed ^ 0xB1A6ULL ^ n);
    auto real = summarize(std::move(real_raw), seed ^ 0x7EA1ULL ^ n);
    const double plan_saved = legacy.median - planned.median;
    const double real_saved = planned.median - real.median;
    const double extra_real_setup = real_setup.median - complex_setup.median;
    return {n, samples, iters, std::move(complex_setup), std::move(real_setup), std::move(legacy), std::move(planned), std::move(real),
            legacy.median / planned.median, planned.median / real.median,
            plan_saved > 0 ? complex_setup.median / plan_saved : std::numeric_limits<double>::infinity(),
            real_saved > 0 ? std::max(0.0, extra_real_setup) / real_saved : std::numeric_limits<double>::infinity()};
}

void planned_tests() {
    std::size_t checks = 0; auto req = [&](bool ok, const char* what) { ++checks; if (!ok) throw std::runtime_error(std::string("planned self-test failure: ") + what); };
    std::mt19937_64 rng(0xC001D00DULL); std::uniform_real_distribution<double> dist(-1.0, 1.0);
    for (std::size_t n : {2u, 4u, 8u, 16u, 64u, 256u}) {
        Radix2Plan plan(n); Vector x(n); for (auto& z : x) z = {dist(rng), dist(rng)}; auto reference = x; radix2_inplace(reference);
        auto y = x; plan.forward_inplace(y); for (std::size_t i = 0; i < n; ++i) req(std::abs(y[i] - reference[i]) < 1e-11 * (1 + std::abs(reference[i])), "planned complex forward");
        plan.inverse_inplace(y); for (std::size_t i = 0; i < n; ++i) req(std::abs(y[i] - x[i]) < 1e-11 * (1 + std::abs(x[i])), "planned complex roundtrip");
        RealVector r(n); for (auto& value : r) value = dist(rng); RealRadix2Plan rp(n); HalfSpectrum half(n / 2 + 1); Vector scratch(n / 2); rp.forward(r, half, scratch);
        Vector full_input(n); for (std::size_t i = 0; i < n; ++i) full_input[i] = {r[i], 0}; radix2_inplace(full_input);
        for (std::size_t k = 0; k <= n / 2; ++k) req(std::abs(half[k] - full_input[k]) < 1e-11 * (1 + std::abs(full_input[k])), "real forward spectrum");
        RealVector restored(n); rp.inverse(half, restored, scratch); for (std::size_t i = 0; i < n; ++i) req(std::abs(restored[i] - r[i]) < 1e-11 * (1 + std::abs(r[i])), "real roundtrip");
    }
    bool threw = false; try { Radix2Plan bad(12); } catch (const std::invalid_argument&) { threw = true; } req(threw, "plan rejects non-power-of-two");
    threw = false; try { RealRadix2Plan bad(1); } catch (const std::invalid_argument&) { threw = true; } req(threw, "real plan rejects N < 2");
    std::cout << "PASS: " << checks << " planned checks\n";
}
}
