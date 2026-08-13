#pragma once

#include "fftlab/fft.hpp"

#include <algorithm>
#include <cstddef>
#include <span>
#include <stdexcept>
#include <vector>

namespace fftlab {

template <FftScalar T = double>
class Radix2Plan {
public:
    explicit Radix2Plan(std::size_t n) : n_(n), bit_reverse_(n), twiddles_(n / 2) {
        if (n == 0 || !pow2(n)) throw std::invalid_argument("Radix2Plan requires power-of-two N >= 1");
        const auto bits = ilog2(n);
        for (std::size_t i = 0; i < n; ++i) {
            auto x = i;
            std::size_t reversed = 0;
            for (std::size_t b = 0; b < bits; ++b) { reversed = (reversed << 1) | (x & 1U); x >>= 1; }
            bit_reverse_[i] = reversed;
        }
        for (std::size_t k = 0; k < n / 2; ++k) {
            const auto angle = -T{2} * std::numbers::pi_v<T> * static_cast<T>(k) / static_cast<T>(n);
            twiddles_[k] = root<T>(angle);
        }
    }

    [[nodiscard]] std::size_t size() const noexcept { return n_; }
    [[nodiscard]] std::size_t scratch_size() const noexcept { return 0; }
    [[nodiscard]] std::size_t stored_twiddles() const noexcept { return twiddles_.size(); }
    [[nodiscard]] std::size_t stored_indices() const noexcept { return bit_reverse_.size(); }

    void forward_inplace(std::span<ComplexT<T>> data) const { execute(data, Direction::Forward); }
    void inverse_inplace(std::span<ComplexT<T>> data) const { execute(data, Direction::Inverse); }

    void forward(const std::span<const ComplexT<T>> input, std::span<ComplexT<T>> output) const {
        copy_checked(input, output); forward_inplace(output);
    }
    void inverse(const std::span<const ComplexT<T>> input, std::span<ComplexT<T>> output) const {
        copy_checked(input, output); inverse_inplace(output);
    }

private:
    void copy_checked(std::span<const ComplexT<T>> input, std::span<ComplexT<T>> output) const {
        if (input.size() != n_ || output.size() != n_) throw std::invalid_argument("Radix2Plan buffer size mismatch");
        std::copy(input.begin(), input.end(), output.begin());
    }

    void execute(std::span<ComplexT<T>> data, Direction direction) const {
        if (data.size() != n_) throw std::invalid_argument("Radix2Plan buffer size mismatch");
        if (n_ <= 1) return;
        for (std::size_t i = 0; i < n_; ++i) if (i < bit_reverse_[i]) std::swap(data[i], data[bit_reverse_[i]]);
        for (std::size_t len = 2; len <= n_;) {
            const auto stride = n_ / len;
            for (std::size_t base = 0; base < n_; base += len) {
                for (std::size_t j = 0; j < len / 2; ++j) {
                    auto w = twiddles_[j * stride];
                    if (direction == Direction::Inverse) w = std::conj(w);
                    const auto u = data[base + j];
                    const auto v = data[base + j + len / 2] * w;
                    data[base + j] = u + v;
                    data[base + j + len / 2] = u - v;
                }
            }
            if (len == n_) break;
            len <<= 1;
        }
        if (direction == Direction::Inverse) {
            const auto scale = T{1} / static_cast<T>(n_);
            for (auto& value : data) value *= scale;
        }
    }

    std::size_t n_{};
    std::vector<std::size_t> bit_reverse_;
    VectorT<T> twiddles_;
};
Radix2Plan(std::size_t) -> Radix2Plan<double>;

template <FftScalar T = double>
class RealRadix2Plan {
public:
    explicit RealRadix2Plan(std::size_t n)
        : n_(n), half_(n / 2), half_plan_(n == 1 ? 1 : n / 2), post_twiddles_(n == 1 ? 1 : n / 2 + 1) {
        if (n == 0 || !pow2(n)) throw std::invalid_argument("RealRadix2Plan requires power-of-two N >= 1");
        if (n == 1) { post_twiddles_[0] = {1, 0}; return; }
        for (std::size_t k = 0; k <= half_; ++k) {
            const auto angle = -T{2} * std::numbers::pi_v<T> * static_cast<T>(k) / static_cast<T>(n_);
            post_twiddles_[k] = root<T>(angle);
        }
    }

    [[nodiscard]] std::size_t size() const noexcept { return n_; }
    [[nodiscard]] std::size_t spectrum_size() const noexcept { return n_ == 1 ? 1 : half_ + 1; }
    [[nodiscard]] std::size_t scratch_size() const noexcept { return n_ == 1 ? 0 : half_; }

    void forward(std::span<const T> input, std::span<ComplexT<T>> output,
                 std::span<ComplexT<T>> scratch) const {
        check(input.size(), output.size(), scratch.size());
        if (n_ == 1) { output[0] = {input[0], 0}; return; }
        for (std::size_t j = 0; j < half_; ++j) scratch[j] = {input[2 * j], input[2 * j + 1]};
        half_plan_.forward_inplace(scratch.first(half_));
        for (std::size_t k = 0; k <= half_; ++k) {
            const auto a = scratch[k % half_];
            const auto b = std::conj(scratch[(half_ - k) % half_]);
            const auto even = T{0.5} * (a + b);
            const auto odd = ComplexT<T>{0, T{-0.5}} * (a - b);
            output[k] = even + post_twiddles_[k] * odd;
        }
    }

    void inverse(std::span<const ComplexT<T>> input, std::span<T> output,
                 std::span<ComplexT<T>> scratch) const {
        check(output.size(), input.size(), scratch.size());
        if (n_ == 1) { output[0] = input[0].real(); return; }
        for (std::size_t k = 0; k < half_; ++k) {
            const auto xk = input[k];
            const auto mirror = std::conj(input[k != 0 ? half_ - k : half_]);
            const auto even = T{0.5} * (xk + mirror);
            const auto odd = T{0.5} * (xk - mirror) / post_twiddles_[k];
            scratch[k] = even + ComplexT<T>{0, 1} * odd;
        }
        half_plan_.inverse_inplace(scratch.first(half_));
        for (std::size_t j = 0; j < half_; ++j) {
            output[2 * j] = scratch[j].real();
            output[2 * j + 1] = scratch[j].imag();
        }
    }

private:
    void check(std::size_t time_size, std::size_t frequency_size, std::size_t scratch_size) const {
        if (time_size != n_ || frequency_size != spectrum_size() || scratch_size < this->scratch_size())
            throw std::invalid_argument("RealRadix2Plan buffer size mismatch");
    }

    std::size_t n_{}, half_{};
    Radix2Plan<T> half_plan_;
    VectorT<T> post_twiddles_;
};
RealRadix2Plan(std::size_t) -> RealRadix2Plan<double>;

} // namespace fftlab
