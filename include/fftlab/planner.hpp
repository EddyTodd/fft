#pragma once

#include "fftlab/convolution_plan.hpp"

namespace fftlab {

template <FftScalar T = double>
class Plan {
    using Storage = std::variant<std::monostate, Radix2Plan<T>, MixedRadixPlan<T>, GoodThomasPlan<T>, RaderPlan<T>,
                                 BluesteinPlan<T>>;

public:
    explicit Plan(std::size_t n, PlanOptions options = {})
        : n_(n), storage_(make_storage(n, options, algorithm_)) {
        scratch_size_ = compute_scratch_size(storage_);
        if (algorithm_ == PlanAlgorithm::Rader || algorithm_ == PlanAlgorithm::Bluestein) {
            inplace_scratch_size_ = scratch_size_ + n_;
        } else if (algorithm_ == PlanAlgorithm::Identity || algorithm_ == PlanAlgorithm::Radix2) {
            inplace_scratch_size_ = 0;
        } else {
            inplace_scratch_size_ = scratch_size_;
        }
    }

    Plan(Plan&&) noexcept = default;
    Plan& operator=(Plan&&) noexcept = default;
    Plan(const Plan&) = delete;
    Plan& operator=(const Plan&) = delete;

    [[nodiscard]] std::size_t size() const noexcept { return n_; }
    [[nodiscard]] PlanAlgorithm algorithm() const noexcept { return algorithm_; }
    [[nodiscard]] std::string_view algorithm_name() const noexcept { return plan_name(algorithm_); }
    [[nodiscard]] std::size_t scratch_size() const noexcept { return scratch_size_; }
    [[nodiscard]] std::size_t inplace_scratch_size() const noexcept { return inplace_scratch_size_; }

    void forward(std::span<const ComplexT<T>> input, std::span<ComplexT<T>> output,
                 std::span<ComplexT<T>> scratch = {}) const {
        execute(input, output, scratch, Direction::Forward);
    }

    void inverse(std::span<const ComplexT<T>> input, std::span<ComplexT<T>> output,
                 std::span<ComplexT<T>> scratch = {}) const {
        execute(input, output, scratch, Direction::Inverse);
    }

    void forward_inplace(std::span<ComplexT<T>> data, std::span<ComplexT<T>> scratch = {}) const {
        execute_inplace(data, scratch, Direction::Forward);
    }

    void inverse_inplace(std::span<ComplexT<T>> data, std::span<ComplexT<T>> scratch = {}) const {
        execute_inplace(data, scratch, Direction::Inverse);
    }

private:
    [[nodiscard]] static std::size_t compute_scratch_size(const Storage& storage) noexcept {
        return std::visit(
            [](const auto& plan) -> std::size_t {
                using P = std::decay_t<decltype(plan)>;
                if constexpr (std::is_same_v<P, std::monostate> || std::is_same_v<P, Radix2Plan<T>>) {
                    return 0;
                } else {
                    return plan.scratch_size();
                }
            },
            storage);
    }

    static Storage make_storage(std::size_t n, PlanOptions options, PlanAlgorithm& algorithm) {
        if (n <= 1) {
            algorithm = PlanAlgorithm::Identity;
            return {};
        }
        const auto caps = plan_capabilities(n);
        const auto force = options.preference;
        if (force == PlanPreference::Radix2) {
            if (!caps.radix2) throw std::invalid_argument("forced radix2 plan unsupported for N");
            algorithm = PlanAlgorithm::Radix2;
            return Storage{std::in_place_type<Radix2Plan<T>>, n};
        }
        if (force == PlanPreference::MixedRadix) {
            if (!caps.mixed_radix) throw std::invalid_argument("forced mixed-radix plan requires composite N");
            algorithm = PlanAlgorithm::MixedRadix;
            return Storage{std::in_place_type<MixedRadixPlan<T>>, n};
        }
        if (force == PlanPreference::GoodThomas) {
            if (!caps.good_thomas) throw std::invalid_argument("forced Good-Thomas plan requires coprime factors");
            algorithm = PlanAlgorithm::GoodThomas;
            return Storage{std::in_place_type<GoodThomasPlan<T>>, n};
        }
        if (force == PlanPreference::Rader) {
            if (!caps.rader) throw std::invalid_argument("forced Rader plan requires prime N >= 3");
            algorithm = PlanAlgorithm::Rader;
            return Storage{std::in_place_type<RaderPlan<T>>, n};
        }
        if (force == PlanPreference::Bluestein) {
            algorithm = PlanAlgorithm::Bluestein;
            return Storage{std::in_place_type<BluesteinPlan<T>>, n};
        }
        if (caps.radix2) {
            algorithm = PlanAlgorithm::Radix2;
            return Storage{std::in_place_type<Radix2Plan<T>>, n};
        }
        if (!is_prime(n)) {
            if (options.prefer_good_thomas && caps.good_thomas) {
                algorithm = PlanAlgorithm::GoodThomas;
                return Storage{std::in_place_type<GoodThomasPlan<T>>, n};
            }
            if (caps.mixed_radix) {
                algorithm = PlanAlgorithm::MixedRadix;
                return Storage{std::in_place_type<MixedRadixPlan<T>>, n};
            }
            algorithm = PlanAlgorithm::Bluestein;
            return Storage{std::in_place_type<BluesteinPlan<T>>, n};
        }
        const auto rader_m = pow2(n - 1) ? n - 1 : next_pow2(2 * (n - 1) - 1);
        const auto blue_m = next_pow2(2 * n - 1);
        if (rader_m < blue_m) {
            algorithm = PlanAlgorithm::Rader;
            return Storage{std::in_place_type<RaderPlan<T>>, n};
        }
        algorithm = PlanAlgorithm::Bluestein;
        return Storage{std::in_place_type<BluesteinPlan<T>>, n};
    }

    void execute(std::span<const ComplexT<T>> input, std::span<ComplexT<T>> output, std::span<ComplexT<T>> scratch,
                 Direction direction) const {
        if (input.size() != n_ || output.size() != n_ || scratch.size() < scratch_size_) {
            throw std::invalid_argument("Plan buffer size mismatch");
        }
        if (n_ <= 1) {
            if (n_ == 1) output[0] = input[0];
            return;
        }
        std::visit(
            [&](const auto& plan) {
                using P = std::decay_t<decltype(plan)>;
                if constexpr (std::is_same_v<P, std::monostate>) {
                    return;
                } else if constexpr (std::is_same_v<P, Radix2Plan<T>>) {
                    if (direction == Direction::Forward) {
                        plan.forward(input, output);
                    } else {
                        plan.inverse(input, output);
                    }
                } else {
                    if (direction == Direction::Forward) {
                        plan.forward(input, output, scratch);
                    } else {
                        plan.inverse(input, output, scratch);
                    }
                }
            },
            storage_);
    }

    void execute_inplace(std::span<ComplexT<T>> data, std::span<ComplexT<T>> scratch, Direction direction) const {
        if (data.size() != n_ || scratch.size() < inplace_scratch_size_) {
            throw std::invalid_argument("Plan buffer size mismatch");
        }
        if (n_ <= 1) return;
        std::visit(
            [&](const auto& plan) {
                using P = std::decay_t<decltype(plan)>;
                if constexpr (std::is_same_v<P, std::monostate>) {
                    return;
                } else if constexpr (std::is_same_v<P, Radix2Plan<T>>) {
                    if (direction == Direction::Forward) {
                        plan.forward_inplace(data);
                    } else {
                        plan.inverse_inplace(data);
                    }
                } else if constexpr (std::is_same_v<P, MixedRadixPlan<T>> || std::is_same_v<P, GoodThomasPlan<T>>) {
                    if (direction == Direction::Forward) {
                        plan.forward_inplace(data, scratch);
                    } else {
                        plan.inverse_inplace(data, scratch);
                    }
                } else {
                    // Rader/Bluestein are naturally out-of-place reductions. Use the
                    // first N scratch elements as a temporary output only when enough
                    // workspace exists; callers wanting minimum scratch should use the
                    // explicit out-of-place API.
                    auto temp = scratch.first(n_);
                    auto algorithm_scratch = scratch.subspan(n_, scratch_size_);
                    if (direction == Direction::Forward) {
                        plan.forward(data, temp, algorithm_scratch);
                    } else {
                        plan.inverse(data, temp, algorithm_scratch);
                    }
                    std::copy(temp.begin(), temp.end(), data.begin());
                }
            },
            storage_);
    }

    std::size_t n_{};
    PlanAlgorithm algorithm_{PlanAlgorithm::Identity};
    Storage storage_;
    std::size_t scratch_size_{};
    std::size_t inplace_scratch_size_{};
};

Plan(std::size_t, PlanOptions = {}) -> Plan<double>;

}  // namespace fftlab
