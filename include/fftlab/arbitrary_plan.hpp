#pragma once

#include "fftlab/planner.hpp"

#include <cstddef>
#include <span>
#include <stdexcept>
#include <vector>

namespace fftlab {

// Compatibility facade for the historical binary64 arbitrary-plan API.
// New code should prefer Plan<T>, which also exposes mixed-radix and PFA plans.
enum class ArbitraryPlanAlgorithm { Radix2, MixedRadix, GoodThomas, Rader, Bluestein };
enum class ArbitraryPlanPolicy { Auto, Bluestein, Rader };

class ArbitraryPlan {
public:
    explicit ArbitraryPlan(std::size_t n, ArbitraryPlanPolicy policy = ArbitraryPlanPolicy::Auto)
        : plan_(n, options_for(policy)) {}

    [[nodiscard]] std::size_t size() const noexcept { return plan_.size(); }
    [[nodiscard]] std::size_t scratch_size() const noexcept { return plan_.scratch_size(); }
    [[nodiscard]] ArbitraryPlanAlgorithm algorithm() const noexcept {
        switch (plan_.algorithm()) {
            case PlanAlgorithm::Radix2: return ArbitraryPlanAlgorithm::Radix2;
            case PlanAlgorithm::MixedRadix: return ArbitraryPlanAlgorithm::MixedRadix;
            case PlanAlgorithm::GoodThomas: return ArbitraryPlanAlgorithm::GoodThomas;
            case PlanAlgorithm::Rader: return ArbitraryPlanAlgorithm::Rader;
            case PlanAlgorithm::Bluestein: return ArbitraryPlanAlgorithm::Bluestein;
            case PlanAlgorithm::Identity: return ArbitraryPlanAlgorithm::Radix2;
        }
        return ArbitraryPlanAlgorithm::Bluestein;
    }

    void forward(const Vector& input, Vector& output, Vector& scratch) const {
        plan_.forward(input, output, scratch);
    }
    void inverse(const Vector& input, Vector& output, Vector& scratch) const {
        plan_.inverse(input, output, scratch);
    }

private:
    [[nodiscard]] static PlanOptions options_for(ArbitraryPlanPolicy policy) {
        switch (policy) {
            case ArbitraryPlanPolicy::Auto: return {};
            case ArbitraryPlanPolicy::Bluestein: return {PlanPreference::Bluestein, true};
            case ArbitraryPlanPolicy::Rader: return {PlanPreference::Rader, true};
        }
        throw std::invalid_argument("unknown arbitrary plan policy");
    }

    Plan<double> plan_;
};

} // namespace fftlab
