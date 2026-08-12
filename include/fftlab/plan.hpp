#pragma once

#include "fftlab/fft.hpp"

#include <cstddef>
#include <cstdint>
#include <string_view>
#include <vector>

namespace fftlab {
using RealVector = std::vector<double>;
using HalfSpectrum = std::vector<Complex>;

class Radix2Plan {
public:
    explicit Radix2Plan(std::size_t n);
    [[nodiscard]] std::size_t size() const noexcept { return n_; }
    [[nodiscard]] std::size_t stored_twiddles() const noexcept { return twiddles_.size(); }
    [[nodiscard]] std::size_t stored_indices() const noexcept { return bit_reverse_.size(); }
    void forward_inplace(Vector& data) const;
    void inverse_inplace(Vector& data) const;

private:
    void execute(Vector& data, bool inverse) const;
    std::size_t n_{};
    std::vector<std::size_t> bit_reverse_;
    Vector twiddles_;
};

class RealRadix2Plan {
public:
    explicit RealRadix2Plan(std::size_t n);
    [[nodiscard]] std::size_t size() const noexcept { return n_; }
    [[nodiscard]] std::size_t spectrum_size() const noexcept { return half_ + 1; }
    [[nodiscard]] std::size_t scratch_size() const noexcept { return half_; }
    void forward(const RealVector& input, HalfSpectrum& output, Vector& scratch) const;
    void inverse(const HalfSpectrum& input, RealVector& output, Vector& scratch) const;

private:
    std::size_t n_{}, half_{};
    Radix2Plan half_plan_;
    Vector post_twiddles_;
};

struct PlanDistribution {
    double min{}, p05{}, median{}, p95{}, max{}, mad{}, ci_lo{}, ci_hi{};
    std::vector<double> raw;
};

struct PlanBenchmark {
    std::size_t n{}, samples{}, iterations_per_sample{};
    PlanDistribution complex_setup, real_setup, legacy_complex, planned_complex, planned_real;
    double plan_speedup{}, real_speedup{}, complex_setup_break_even_transforms{}, real_setup_break_even_transforms{};
};

PlanBenchmark benchmark_plans(std::size_t n, std::size_t samples = 31, std::size_t warmups = 5,
                              double target_ms = 5.0, std::uint64_t seed = 0xF17F17ULL);
void planned_tests();
}
