#include "fftlab/oracle.hpp"
#include "fftlab/planner.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <numbers>
#include <random>
#include <stdexcept>
#include <string>
#include <string_view>

using namespace fftlab;

namespace {

std::size_t checks = 0;

void require(bool condition, std::string_view message) {
    ++checks;
    if (!condition) {
        throw std::runtime_error(std::string(message));
    }
}

template <FftScalar T>
[[nodiscard]] constexpr T tolerance() {
    return std::same_as<T, float> ? static_cast<T>(8e-4) : static_cast<T>(8e-10);
}

template <FftScalar T>
[[nodiscard]] bool approx(ComplexT<T> actual, ComplexT<T> reference, T multiplier = T{1}) {
    return std::abs(actual - reference) <=
           tolerance<T>() * multiplier * (T{1} + std::abs(reference));
}

template <FftScalar T>
[[nodiscard]] bool approx_real(T actual, T reference, T multiplier = T{1}) {
    return std::abs(actual - reference) <=
           tolerance<T>() * multiplier * (T{1} + std::abs(reference));
}

template <FftScalar T>
void require_oracle_close(ComplexT<T> actual, OracleComplex reference, T multiplier,
                          std::string_view message) {
    const OracleComplex converted{static_cast<long double>(actual.real()),
                                  static_cast<long double>(actual.imag())};
    const auto limit = static_cast<long double>(tolerance<T>() * multiplier) *
                       (1.0L + std::abs(reference));
    require(std::abs(converted - reference) <= limit, message);
}

template <FftScalar T>
void oracle_matrix() {
    std::mt19937_64 rng{0xF17F17};
    std::uniform_real_distribution<double> distribution{-1.0, 1.0};

    for (std::size_t n : {0u, 1u, 2u, 3u, 4u, 5u, 6u, 7u, 8u, 9u, 10u, 12u, 15u,
                          16u, 17u, 18u, 20u, 21u, 25u, 27u, 31u, 35u, 49u, 60u, 64u,
                          75u, 121u, 127u}) {
        VectorT<T> input(n);
        for (auto& value : input) {
            value = {static_cast<T>(distribution(rng)), static_cast<T>(distribution(rng))};
        }

        const auto reference = oracle_dft<T>(input);
        const auto automatic = transform(input);
        for (std::size_t i = 0; i < n; ++i) {
            require_oracle_close(automatic[i], reference[i], T{4}, "auto algorithm oracle");
        }

        Plan<T> plan(n);
        VectorT<T> output(n);
        VectorT<T> scratch(plan.scratch_size());
        plan.forward(input, output, scratch);
        for (std::size_t i = 0; i < n; ++i) {
            require_oracle_close(output[i], reference[i], T{4}, "structural plan oracle");
        }

        VectorT<T> roundtrip(n);
        plan.inverse(output, roundtrip, scratch);
        for (std::size_t i = 0; i < n; ++i) {
            require(approx(roundtrip[i], input[i], T{3}), "plan roundtrip");
        }
    }
}

template <FftScalar T>
void catalog_oracle() {
    std::mt19937_64 rng{0xCA7A106};
    std::uniform_real_distribution<double> distribution{-1.0, 1.0};

    for (std::size_t n : {0u, 1u, 2u, 3u, 4u, 5u, 6u, 7u, 8u, 10u, 12u, 15u, 16u,
                          17u, 20u, 21u, 25u, 31u, 32u}) {
        VectorT<T> input(n);
        for (auto& value : input) {
            value = {static_cast<T>(distribution(rng)), static_cast<T>(distribution(rng))};
        }
        const auto reference = oracle_dft<T>(input);

        for (const auto algorithm : all_algorithms) {
            if (algorithm == Algorithm::Auto || algorithm == Algorithm::Dft ||
                !supports_algorithm(algorithm, n)) {
                continue;
            }

            const auto output = transform(input, algorithm, Direction::Forward);
            require(output.size() == n, "catalog output size");
            for (std::size_t i = 0; i < n; ++i) {
                require_oracle_close(output[i], reference[i], T{10}, "catalog forward oracle");
            }

            const auto roundtrip = transform(output, algorithm, Direction::Inverse);
            for (std::size_t i = 0; i < n; ++i) {
                require(approx(roundtrip[i], input[i], T{8}), "catalog inverse roundtrip");
            }
        }
    }
}

template <FftScalar T>
void codelet_tests() {
    for (const std::size_t radix : {2u, 3u, 4u, 5u, 7u}) {
        SmallDftCodelet<T> codelet(radix);
        require(codelet.radix() == radix, "codelet radix");
        if (radix == 2 || radix == 4) {
            require(codelet.stored_roots() == 0, "specialized codelet has no root matrix");
        } else {
            require(codelet.stored_roots() == radix * radix, "planned codelet root matrix");
        }

        VectorT<T> input(radix);
        for (std::size_t i = 0; i < radix; ++i) {
            input[i] = {static_cast<T>(i + 1) / static_cast<T>(radix),
                        static_cast<T>(2 * i + 1) / static_cast<T>(radix + 1)};
        }

        const auto reference = oracle_dft<T>(input);
        auto output = input;
        codelet.execute(output, Direction::Forward);
        for (std::size_t i = 0; i < radix; ++i) {
            require_oracle_close(output[i], reference[i], T{4}, "codelet oracle");
        }

        codelet.execute(output, Direction::Inverse);
        for (std::size_t i = 0; i < radix; ++i) {
            require(approx(output[i], input[i], T{4}), "codelet roundtrip");
        }
    }
}

template <FftScalar T>
void identities() {
    for (std::size_t n : {3u, 8u, 12u, 15u, 17u, 25u, 31u, 64u}) {
        VectorT<T> impulse(n);
        impulse[0] = {1, 0};
        auto output = transform(impulse);
        for (const auto value : output) {
            require(approx(value, {1, 0}, T{2}), "impulse identity");
        }

        VectorT<T> constant(n, {T{1}, T{0}});
        output = transform(constant);
        require(approx(output[0], {static_cast<T>(n), 0}, T{2}), "constant DC");
        for (std::size_t k = 1; k < n; ++k) {
            require(std::abs(output[k]) <= tolerance<T>() * T{8} * static_cast<T>(n),
                    "constant non-DC");
        }

        const std::size_t tone = n > 1 ? std::min<std::size_t>(2, n - 1) : 0;
        VectorT<T> single_tone(n);
        for (std::size_t j = 0; j < n; ++j) {
            const T angle = T{2} * std::numbers::pi_v<T> * static_cast<T>(tone * j) /
                            static_cast<T>(n);
            single_tone[j] = root<T>(angle);
        }
        output = transform(single_tone);
        for (std::size_t k = 0; k < n; ++k) {
            if (k == tone) {
                require(approx(output[k], {static_cast<T>(n), 0}, T{3}), "single tone bin");
            } else {
                require(std::abs(output[k]) <= tolerance<T>() * T{12} * static_cast<T>(n),
                        "single tone leakage");
            }
        }
    }
}

template <FftScalar T>
void real_tests() {
    std::mt19937_64 rng{2026};
    std::uniform_real_distribution<double> distribution{-1.0, 1.0};

    for (std::size_t n : {1u, 2u, 4u, 8u, 16u, 64u}) {
        RealRadix2Plan<T> plan(n);
        std::vector<T> input(n);
        for (auto& value : input) {
            value = static_cast<T>(distribution(rng));
        }

        std::vector<ComplexT<T>> spectrum(plan.spectrum_size());
        std::vector<ComplexT<T>> scratch(plan.scratch_size());
        plan.forward(input, spectrum, scratch);

        VectorT<T> complex_input(n);
        for (std::size_t i = 0; i < n; ++i) {
            complex_input[i] = {input[i], 0};
        }
        const auto full = dft(complex_input);
        for (std::size_t k = 0; k < spectrum.size(); ++k) {
            require(approx(spectrum[k], full[k], T{3}), "real half spectrum");
        }
        if (n > 1) {
            for (std::size_t k = 1; k < n / 2; ++k) {
                require(approx(full[n - k], std::conj(full[k]), T{3}), "Hermitian symmetry");
            }
        }

        std::vector<T> roundtrip(n);
        plan.inverse(spectrum, roundtrip, scratch);
        for (std::size_t i = 0; i < n; ++i) {
            require(approx_real(roundtrip[i], input[i], T{3}), "real roundtrip");
        }
    }
}

template <FftScalar T>
void explicit_mechanisms() {
    std::mt19937_64 rng{99};
    std::uniform_real_distribution<double> distribution{-1.0, 1.0};

    for (std::size_t n : {12u, 18u, 20u, 25u, 27u, 45u, 49u, 60u, 75u, 121u}) {
        VectorT<T> input(n);
        for (auto& value : input) {
            value = {static_cast<T>(distribution(rng)), static_cast<T>(distribution(rng))};
        }
        const auto reference = oracle_dft<T>(input);

        MixedRadixPlan<T> plan(n);
        auto output = input;
        VectorT<T> scratch(plan.scratch_size());
        plan.forward_inplace(output, scratch);
        for (std::size_t i = 0; i < n; ++i) {
            require_oracle_close(output[i], reference[i], T{6}, "planned mixed radix");
        }
    }

    for (std::size_t n : {6u, 10u, 12u, 15u, 18u, 20u, 21u, 35u, 45u, 60u, 75u, 143u}) {
        const auto split = coprime_factor_split(n);
        if (split.first == 0) {
            continue;
        }

        VectorT<T> input(n);
        for (auto& value : input) {
            value = {static_cast<T>(distribution(rng)), static_cast<T>(distribution(rng))};
        }
        const auto reference = oracle_dft<T>(input);

        GoodThomasPlan<T> plan(n);
        require(plan.twiddle_count() == 0, "PFA top-level is twiddle free");
        auto output = input;
        VectorT<T> scratch(plan.scratch_size());
        plan.forward_inplace(output, scratch);
        for (std::size_t i = 0; i < n; ++i) {
            require_oracle_close(output[i], reference[i], T{8}, "Good Thomas plan");
        }
    }

    for (std::size_t n : {8u, 16u, 32u, 64u, 128u}) {
        VectorT<T> input(n);
        for (auto& value : input) {
            value = {static_cast<T>(distribution(rng)), static_cast<T>(distribution(rng))};
        }
        const auto reference = oracle_dft<T>(input);
        const auto output = modified_split_radix(input);
        for (std::size_t i = 0; i < n; ++i) {
            require_oracle_close(output[i], reference[i], T{8}, "modified split radix");
        }
    }
}

void planner_policy() {
    require(Plan<double>(0).algorithm() == PlanAlgorithm::Identity, "N0 identity");
    require(Plan<double>(1).algorithm() == PlanAlgorithm::Identity, "N1 identity");
    require(Plan<double>(64).algorithm() == PlanAlgorithm::Radix2, "power2 radix2");
    require(Plan<double>(25).algorithm() == PlanAlgorithm::MixedRadix, "prime power mixed");
    require(Plan<double>(15).algorithm() == PlanAlgorithm::GoodThomas, "coprime PFA");
    require(Plan<double>(17).algorithm() == PlanAlgorithm::Rader, "short Rader");
    require(Plan<double>(31).algorithm() == PlanAlgorithm::Bluestein,
            "equal-convolution prime Bluestein");
    require(Plan<double>(121).algorithm() == PlanAlgorithm::Bluestein,
            "rough composite Bluestein fallback");

    bool threw = false;
    try {
        Plan<double> invalid(12, {PlanPreference::Rader, true});
        (void)invalid;
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    require(threw, "forced invalid policy rejects");
}

}  // namespace

int main() {
    oracle_matrix<float>();
    oracle_matrix<double>();
    catalog_oracle<float>();
    catalog_oracle<double>();
    codelet_tests<float>();
    codelet_tests<double>();
    identities<float>();
    identities<double>();
    real_tests<float>();
    real_tests<double>();
    explicit_mechanisms<float>();
    explicit_mechanisms<double>();
    planner_policy();
    std::cout << "PASS: " << checks << " checks\n";
}
