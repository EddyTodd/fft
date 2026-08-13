#pragma once

#include "fftlab/mixed_plan.hpp"

namespace fftlab {

template <FftScalar T = double>
class GoodThomasPlan {
public:
    explicit GoodThomasPlan(std::size_t n) : GoodThomasPlan(n, coprime_factor_split(n)) {}
    GoodThomasPlan(std::size_t n, std::size_t a, std::size_t b) : GoodThomasPlan(n, std::pair{a, b}) {}

    [[nodiscard]] std::size_t size() const noexcept { return n_; }
    [[nodiscard]] std::size_t factor_a() const noexcept { return a_; }
    [[nodiscard]] std::size_t factor_b() const noexcept { return b_; }
    [[nodiscard]] std::size_t scratch_size() const noexcept { return n_ + 2 * std::max(a_, b_); }
    [[nodiscard]] std::size_t twiddle_count() const noexcept { return 0; }

    void forward_inplace(std::span<ComplexT<T>> data, std::span<ComplexT<T>> scratch) const {
        execute(data, scratch, Direction::Forward);
    }
    void inverse_inplace(std::span<ComplexT<T>> data, std::span<ComplexT<T>> scratch) const {
        execute(data, scratch, Direction::Inverse);
    }
    void forward(std::span<const ComplexT<T>> input, std::span<ComplexT<T>> output,
                 std::span<ComplexT<T>> scratch) const {
        copy_checked(input, output, scratch); forward_inplace(output, scratch);
    }
    void inverse(std::span<const ComplexT<T>> input, std::span<ComplexT<T>> output,
                 std::span<ComplexT<T>> scratch) const {
        copy_checked(input, output, scratch); inverse_inplace(output, scratch);
    }

private:
    GoodThomasPlan(std::size_t n, std::pair<std::size_t, std::size_t> factors)
        : n_(n), a_(factors.first), b_(factors.second), a_plan_(a_), b_plan_(b_), input_map_(n), output_map_(n) {
        if (a_ <= 1 || b_ <= 1 || a_ * b_ != n_ || std::gcd(a_, b_) != 1)
            throw std::invalid_argument("GoodThomasPlan requires coprime N=a*b factors");
        const auto inv_b_mod_a = modular_inverse(b_ % a_, a_);
        const auto inv_a_mod_b = modular_inverse(a_ % b_, b_);
        for (std::size_t n1 = 0; n1 < a_; ++n1) {
            for (std::size_t n2 = 0; n2 < b_; ++n2) {
                const auto matrix_index = n1 * b_ + n2;
                input_map_[matrix_index] = (mul_mod(mul_mod(n1, b_, n_), inv_b_mod_a, n_) +
                                            mul_mod(mul_mod(n2, a_, n_), inv_a_mod_b, n_)) % n_;
            }
        }
        for (std::size_t k1 = 0; k1 < a_; ++k1)
            for (std::size_t k2 = 0; k2 < b_; ++k2)
                output_map_[k1 * b_ + k2] = (k1 * b_ + k2 * a_) % n_;
    }

    void execute(std::span<ComplexT<T>> data, std::span<ComplexT<T>> scratch, Direction direction) const {
        if (data.size() != n_ || scratch.size() < scratch_size()) throw std::invalid_argument("GoodThomasPlan buffer size mismatch");
        auto matrix = scratch.first(n_);
        auto buffer = scratch.subspan(n_, 2 * std::max(a_, b_));
        auto work = buffer.first(std::max(a_, b_));
        auto work_scratch = buffer.subspan(std::max(a_, b_), std::max(a_, b_));
        for (std::size_t i = 0; i < n_; ++i) matrix[i] = data[input_map_[i]];

        for (std::size_t n1 = 0; n1 < a_; ++n1) {
            auto row = matrix.subspan(n1 * b_, b_);
            if (direction == Direction::Forward) b_plan_.forward_inplace(row, work_scratch.first(b_));
            else b_plan_.inverse_inplace(row, work_scratch.first(b_));
        }
        for (std::size_t k2 = 0; k2 < b_; ++k2) {
            for (std::size_t n1 = 0; n1 < a_; ++n1) work[n1] = matrix[n1 * b_ + k2];
            if (direction == Direction::Forward) a_plan_.forward_inplace(work.first(a_), work_scratch.first(a_));
            else a_plan_.inverse_inplace(work.first(a_), work_scratch.first(a_));
            for (std::size_t k1 = 0; k1 < a_; ++k1) matrix[k1 * b_ + k2] = work[k1];
        }
        for (std::size_t i = 0; i < n_; ++i) data[output_map_[i]] = matrix[i];
    }

    void copy_checked(std::span<const ComplexT<T>> input, std::span<ComplexT<T>> output,
                      std::span<ComplexT<T>> scratch) const {
        if (input.size() != n_ || output.size() != n_ || scratch.size() < scratch_size())
            throw std::invalid_argument("GoodThomasPlan buffer size mismatch");
        std::copy(input.begin(), input.end(), output.begin());
    }

    std::size_t n_{}, a_{}, b_{};
    MixedRadixPlan<T> a_plan_;
    MixedRadixPlan<T> b_plan_;
    std::vector<std::size_t> input_map_, output_map_;
};
GoodThomasPlan(std::size_t) -> GoodThomasPlan<double>;
GoodThomasPlan(std::size_t, std::size_t, std::size_t) -> GoodThomasPlan<double>;

} // namespace fftlab
