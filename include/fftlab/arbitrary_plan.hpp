#pragma once

#include "fftlab/plan.hpp"

#include <cstddef>
#include <variant>
#include <vector>

namespace fftlab {

enum class ArbitraryPlanAlgorithm { Radix2, Bluestein, Rader };
enum class ArbitraryPlanPolicy { Auto, Bluestein, Rader };

class BluesteinPlan {
public:
    explicit BluesteinPlan(std::size_t n);

    [[nodiscard]] std::size_t size() const noexcept { return n_; }
    [[nodiscard]] std::size_t convolution_size() const noexcept { return m_; }
    [[nodiscard]] std::size_t scratch_size() const noexcept { return m_; }
    [[nodiscard]] std::size_t persistent_complex_values() const noexcept {
        return chirp_.size() + kernel_spectrum_.size() + convolution_plan_.stored_twiddles();
    }
    [[nodiscard]] std::size_t persistent_indices() const noexcept {
        return convolution_plan_.stored_indices();
    }

    void forward(const Vector& input, Vector& output, Vector& scratch) const;
    void inverse(const Vector& input, Vector& output, Vector& scratch) const;

private:
    void execute(const Vector& input, Vector& output, Vector& scratch, bool inverse) const;

    std::size_t n_{}, m_{};
    Radix2Plan convolution_plan_;
    Vector chirp_;
    Vector kernel_spectrum_;
};

class RaderPlan {
public:
    explicit RaderPlan(std::size_t prime_n);

    [[nodiscard]] std::size_t size() const noexcept { return n_; }
    [[nodiscard]] std::size_t cyclic_size() const noexcept { return l_; }
    [[nodiscard]] std::size_t convolution_size() const noexcept { return m_; }
    [[nodiscard]] std::size_t scratch_size() const noexcept { return m_; }
    [[nodiscard]] bool direct_cyclic_fft() const noexcept { return direct_cyclic_; }
    [[nodiscard]] std::size_t persistent_complex_values() const noexcept {
        return kernel_spectrum_.size() + convolution_plan_.stored_twiddles();
    }
    [[nodiscard]] std::size_t persistent_indices() const noexcept {
        return output_permutation_.size() + input_permutation_.size() + convolution_plan_.stored_indices();
    }

    void forward(const Vector& input, Vector& output, Vector& scratch) const;
    void inverse(const Vector& input, Vector& output, Vector& scratch) const;

private:
    void execute(const Vector& input, Vector& output, Vector& scratch, bool inverse) const;

    std::size_t n_{}, l_{}, m_{};
    bool direct_cyclic_{};
    Radix2Plan convolution_plan_;
    std::vector<std::size_t> output_permutation_;
    std::vector<std::size_t> input_permutation_;
    Vector kernel_spectrum_;
};

class ArbitraryPlan {
public:
    explicit ArbitraryPlan(std::size_t n, ArbitraryPlanPolicy policy = ArbitraryPlanPolicy::Auto);

    [[nodiscard]] std::size_t size() const noexcept { return n_; }
    [[nodiscard]] ArbitraryPlanAlgorithm algorithm() const noexcept { return algorithm_; }
    [[nodiscard]] std::size_t scratch_size() const noexcept;

    void forward(const Vector& input, Vector& output, Vector& scratch) const;
    void inverse(const Vector& input, Vector& output, Vector& scratch) const;

private:
    using Storage = std::variant<Radix2Plan, BluesteinPlan, RaderPlan>;
    static Storage make_storage(std::size_t n, ArbitraryPlanPolicy policy, ArbitraryPlanAlgorithm& algorithm);

    std::size_t n_{};
    ArbitraryPlanAlgorithm algorithm_{ArbitraryPlanAlgorithm::Bluestein};
    Storage storage_;
};

void arbitrary_plan_tests();

} // namespace fftlab
