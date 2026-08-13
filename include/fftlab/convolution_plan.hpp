#pragma once

#include "fftlab/good_thomas_plan.hpp"

namespace fftlab {

template <FftScalar T = double>
class BluesteinPlan {
public:
    explicit BluesteinPlan(std::size_t n)
        : n_(n), m_(workspace(n)), convolution_plan_(m_), chirp_(n), kernel_spectrum_(m_) {
        for (std::size_t k = 0; k < n_; ++k) {
            const auto c = forward_chirp(k, n_);
            chirp_[k] = c;
            const auto b = std::conj(c);
            kernel_spectrum_[k] = b;
            if (k != 0) kernel_spectrum_[m_ - k] = b;
        }
        convolution_plan_.forward_inplace(kernel_spectrum_);
    }
    [[nodiscard]] std::size_t size() const noexcept { return n_; }
    [[nodiscard]] std::size_t convolution_size() const noexcept { return m_; }
    [[nodiscard]] std::size_t scratch_size() const noexcept { return m_; }
    void forward(std::span<const ComplexT<T>> input, std::span<ComplexT<T>> output, std::span<ComplexT<T>> scratch) const {
        execute(input, output, scratch, Direction::Forward);
    }
    void inverse(std::span<const ComplexT<T>> input, std::span<ComplexT<T>> output, std::span<ComplexT<T>> scratch) const {
        execute(input, output, scratch, Direction::Inverse);
    }
private:
    [[nodiscard]] static std::size_t workspace(std::size_t n) {
        if (n == 0) throw std::invalid_argument("BluesteinPlan requires N >= 1");
        if (n == 1) return 1;
        if (n > (std::numeric_limits<std::size_t>::max() / 2) + 1) throw std::length_error("BluesteinPlan workspace overflow");
        return next_pow2(2 * n - 1);
    }
    [[nodiscard]] static ComplexT<T> forward_chirp(std::size_t k, std::size_t n) {
        const long double kd = static_cast<long double>(k);
        const long double phase = std::fmod(kd * kd, 2.0L * static_cast<long double>(n)) / static_cast<long double>(n);
        return root<T>(-std::numbers::pi_v<T> * static_cast<T>(phase));
    }
    void execute(std::span<const ComplexT<T>> input, std::span<ComplexT<T>> output,
                 std::span<ComplexT<T>> scratch, Direction direction) const {
        if (input.size() != n_ || output.size() != n_ || scratch.size() < m_)
            throw std::invalid_argument("BluesteinPlan buffer size mismatch");
        std::fill(scratch.begin(), scratch.begin() + static_cast<std::ptrdiff_t>(m_), ComplexT<T>{});
        for (std::size_t k = 0; k < n_; ++k) {
            const auto x = direction == Direction::Inverse ? std::conj(input[k]) : input[k];
            scratch[k] = x * chirp_[k];
        }
        auto work = scratch.first(m_);
        convolution_plan_.forward_inplace(work);
        for (std::size_t k = 0; k < m_; ++k) work[k] *= kernel_spectrum_[k];
        convolution_plan_.inverse_inplace(work);
        if (direction == Direction::Forward) {
            for (std::size_t k = 0; k < n_; ++k) output[k] = work[k] * chirp_[k];
        } else {
            const auto scale = T{1} / static_cast<T>(n_);
            for (std::size_t k = 0; k < n_; ++k) output[k] = std::conj(work[k] * chirp_[k]) * scale;
        }
    }
    std::size_t n_{}, m_{};
    Radix2Plan<T> convolution_plan_;
    VectorT<T> chirp_, kernel_spectrum_;
};
BluesteinPlan(std::size_t) -> BluesteinPlan<double>;

template <FftScalar T = double>
class RaderPlan {
public:
    explicit RaderPlan(std::size_t n)
        : n_(n), l_(n > 0 ? n - 1 : 0), m_(workspace_checked(n)), direct_cyclic_(m_ == l_),
          convolution_plan_(m_), output_permutation_(l_), input_permutation_(l_), kernel_spectrum_(m_) {
        const auto g = primitive_root_prime(n_);
        output_permutation_[0] = 1;
        for (std::size_t q = 1; q < l_; ++q) output_permutation_[q] = mul_mod(output_permutation_[q - 1], g, n_);
        for (std::size_t q = 0; q < l_; ++q) input_permutation_[q] = output_permutation_[(l_ - q) % l_];
        for (std::size_t q = 0; q < l_; ++q)
            kernel_spectrum_[q] = root<T>(-T{2} * std::numbers::pi_v<T> * static_cast<T>(output_permutation_[q]) / static_cast<T>(n_));
        convolution_plan_.forward_inplace(kernel_spectrum_);
    }
    [[nodiscard]] std::size_t size() const noexcept { return n_; }
    [[nodiscard]] std::size_t cyclic_size() const noexcept { return l_; }
    [[nodiscard]] std::size_t convolution_size() const noexcept { return m_; }
    [[nodiscard]] bool direct_cyclic_fft() const noexcept { return direct_cyclic_; }
    [[nodiscard]] std::size_t scratch_size() const noexcept { return m_; }
    void forward(std::span<const ComplexT<T>> input, std::span<ComplexT<T>> output, std::span<ComplexT<T>> scratch) const {
        execute(input, output, scratch, Direction::Forward);
    }
    void inverse(std::span<const ComplexT<T>> input, std::span<ComplexT<T>> output, std::span<ComplexT<T>> scratch) const {
        execute(input, output, scratch, Direction::Inverse);
    }
private:
    [[nodiscard]] static std::size_t workspace_checked(std::size_t n) {
        if (n < 3 || !is_prime(n)) throw std::invalid_argument("RaderPlan requires prime N >= 3");
        const auto l = n - 1;
        if (pow2(l)) return l;
        if (l > (std::numeric_limits<std::size_t>::max() / 2) + 1) throw std::length_error("RaderPlan workspace overflow");
        return next_pow2(2 * l - 1);
    }
    void execute(std::span<const ComplexT<T>> input, std::span<ComplexT<T>> output,
                 std::span<ComplexT<T>> scratch, Direction direction) const {
        if (input.size() != n_ || output.size() != n_ || scratch.size() < m_)
            throw std::invalid_argument("RaderPlan buffer size mismatch");
        auto work = scratch.first(m_);
        std::fill(work.begin(), work.end(), ComplexT<T>{});
        ComplexT<T> dc{};
        for (std::size_t i = 0; i < n_; ++i) dc += direction == Direction::Inverse ? std::conj(input[i]) : input[i];
        for (std::size_t q = 0; q < l_; ++q)
            work[q] = direction == Direction::Inverse ? std::conj(input[input_permutation_[q]]) : input[input_permutation_[q]];
        convolution_plan_.forward_inplace(work);
        for (std::size_t i = 0; i < m_; ++i) work[i] *= kernel_spectrum_[i];
        convolution_plan_.inverse_inplace(work);
        const auto x0 = direction == Direction::Inverse ? std::conj(input[0]) : input[0];
        if (direction == Direction::Forward) {
            output[0] = dc;
            for (std::size_t q = 0; q < l_; ++q) {
                auto c = work[q];
                if (!direct_cyclic_ && q + l_ < 2 * l_ - 1) c += work[q + l_];
                output[output_permutation_[q]] = x0 + c;
            }
        } else {
            const auto scale = T{1} / static_cast<T>(n_);
            output[0] = std::conj(dc) * scale;
            for (std::size_t q = 0; q < l_; ++q) {
                auto c = work[q];
                if (!direct_cyclic_ && q + l_ < 2 * l_ - 1) c += work[q + l_];
                output[output_permutation_[q]] = std::conj(x0 + c) * scale;
            }
        }
    }
    std::size_t n_{}, l_{}, m_{};
    bool direct_cyclic_{};
    Radix2Plan<T> convolution_plan_;
    std::vector<std::size_t> output_permutation_, input_permutation_;
    VectorT<T> kernel_spectrum_;
};
RaderPlan(std::size_t) -> RaderPlan<double>;

} // namespace fftlab
