#include "fftlab/kernel.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <iostream>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <vector>

#if (defined(__GNUC__) || defined(__clang__)) && (defined(__x86_64__) || defined(__i386__))
#include <immintrin.h>
#define FFTLAB_X86_TARGETS 1
#else
#define FFTLAB_X86_TARGETS 0
#endif

namespace fftlab {
namespace {

using Swap = std::pair<std::size_t, std::size_t>;
using KernelFn = void (*)(Complex*, std::size_t, const Swap*, std::size_t,
                          const std::size_t*, std::size_t, const Complex*, bool);
using Clock = std::chrono::steady_clock;

constexpr double pi = std::numbers::pi_v<double>;

void apply_permutation(Complex* data, const Swap* swaps, std::size_t count) {
    for (std::size_t i = 0; i < count; ++i) {
        std::swap(data[swaps[i].first], data[swaps[i].second]);
    }
}

void scalar_kernel(Complex* data, std::size_t n, const Swap* swaps, std::size_t swap_count,
                   const std::size_t* offsets, std::size_t stage_count,
                   const Complex* twiddles, bool inverse) {
    apply_permutation(data, swaps, swap_count);

    std::size_t stage = 0;
    for (std::size_t len = 2; stage < stage_count; ++stage, len <<= 1) {
        const auto* stage_twiddles = twiddles + offsets[stage];
        const auto half = len / 2;
        for (std::size_t base = 0; base < n; base += len) {
            for (std::size_t j = 0; j < half; ++j) {
                const Complex w = inverse ? std::conj(stage_twiddles[j]) : stage_twiddles[j];
                const Complex u = data[base + j];
                const Complex v = data[base + j + half] * w;
                data[base + j] = u + v;
                data[base + j + half] = u - v;
            }
        }
    }

    if (inverse) {
        const double scale = 1.0 / static_cast<double>(n);
        for (std::size_t i = 0; i < n; ++i) data[i] *= scale;
    }
}

#if FFTLAB_X86_TARGETS
__attribute__((target("avx2,fma")))
void avx2_kernel(Complex* data, std::size_t n, const Swap* swaps, std::size_t swap_count,
                 const std::size_t* offsets, std::size_t stage_count,
                 const Complex* twiddles, bool inverse) {
    apply_permutation(data, swaps, swap_count);
    const __m256d conjugate_mask = _mm256_setr_pd(1.0, -1.0, 1.0, -1.0);

    std::size_t stage = 0;
    for (std::size_t len = 2; stage < stage_count; ++stage, len <<= 1) {
        const auto* stage_twiddles = twiddles + offsets[stage];
        const auto half = len / 2;
        for (std::size_t base = 0; base < n; base += len) {
            std::size_t j = 0;
            for (; j + 2 <= half; j += 2) {
                const auto u = _mm256_loadu_pd(reinterpret_cast<double*>(data + base + j));
                const auto v = _mm256_loadu_pd(reinterpret_cast<double*>(data + base + j + half));
                auto w = _mm256_loadu_pd(reinterpret_cast<const double*>(stage_twiddles + j));
                if (inverse) w = _mm256_mul_pd(w, conjugate_mask);

                const auto wr = _mm256_movedup_pd(w);
                const auto wi = _mm256_permute_pd(w, 0xF);
                const auto swapped_v = _mm256_permute_pd(v, 0x5);
                const auto cross = _mm256_mul_pd(swapped_v, wi);
                const auto product = _mm256_fmaddsub_pd(v, wr, cross);

                _mm256_storeu_pd(reinterpret_cast<double*>(data + base + j),
                                 _mm256_add_pd(u, product));
                _mm256_storeu_pd(reinterpret_cast<double*>(data + base + j + half),
                                 _mm256_sub_pd(u, product));
            }
            for (; j < half; ++j) {
                const Complex w = inverse ? std::conj(stage_twiddles[j]) : stage_twiddles[j];
                const Complex u = data[base + j];
                const Complex v = data[base + j + half] * w;
                data[base + j] = u + v;
                data[base + j + half] = u - v;
            }
        }
    }

    if (inverse) {
        const auto scale = _mm256_set1_pd(1.0 / static_cast<double>(n));
        std::size_t i = 0;
        for (; i + 2 <= n; i += 2) {
            const auto x = _mm256_loadu_pd(reinterpret_cast<double*>(data + i));
            _mm256_storeu_pd(reinterpret_cast<double*>(data + i), _mm256_mul_pd(x, scale));
        }
        for (; i < n; ++i) data[i] /= static_cast<double>(n);
    }
}

__attribute__((target("avx512f,avx512dq,avx512vl,fma")))
void avx512_kernel(Complex* data, std::size_t n, const Swap* swaps, std::size_t swap_count,
                   const std::size_t* offsets, std::size_t stage_count,
                   const Complex* twiddles, bool inverse) {
    apply_permutation(data, swaps, swap_count);
    const __m512d multiply_sign = _mm512_setr_pd(-1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0);
    const __m512d conjugate_mask = _mm512_setr_pd(1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0);

    std::size_t stage = 0;
    for (std::size_t len = 2; stage < stage_count; ++stage, len <<= 1) {
        const auto* stage_twiddles = twiddles + offsets[stage];
        const auto half = len / 2;
        for (std::size_t base = 0; base < n; base += len) {
            std::size_t j = 0;
            for (; j + 4 <= half; j += 4) {
                const auto u = _mm512_loadu_pd(reinterpret_cast<double*>(data + base + j));
                const auto v = _mm512_loadu_pd(reinterpret_cast<double*>(data + base + j + half));
                auto w = _mm512_loadu_pd(reinterpret_cast<const double*>(stage_twiddles + j));
                if (inverse) w = _mm512_mul_pd(w, conjugate_mask);

                const auto wr = _mm512_movedup_pd(w);
                const auto wi = _mm512_permute_pd(w, 0xFF);
                const auto swapped_v = _mm512_permute_pd(v, 0x55);
                const auto cross = _mm512_mul_pd(_mm512_mul_pd(swapped_v, wi), multiply_sign);
                const auto product = _mm512_fmadd_pd(v, wr, cross);

                _mm512_storeu_pd(reinterpret_cast<double*>(data + base + j),
                                 _mm512_add_pd(u, product));
                _mm512_storeu_pd(reinterpret_cast<double*>(data + base + j + half),
                                 _mm512_sub_pd(u, product));
            }
            for (; j < half; ++j) {
                const Complex w = inverse ? std::conj(stage_twiddles[j]) : stage_twiddles[j];
                const Complex u = data[base + j];
                const Complex v = data[base + j + half] * w;
                data[base + j] = u + v;
                data[base + j + half] = u - v;
            }
        }
    }

    if (inverse) {
        const auto scale = _mm512_set1_pd(1.0 / static_cast<double>(n));
        std::size_t i = 0;
        for (; i + 4 <= n; i += 4) {
            const auto x = _mm512_loadu_pd(reinterpret_cast<double*>(data + i));
            _mm512_storeu_pd(reinterpret_cast<double*>(data + i), _mm512_mul_pd(x, scale));
        }
        for (; i < n; ++i) data[i] /= static_cast<double>(n);
    }
}
#endif

KernelFn function_for(KernelIsa isa) {
    if (isa == KernelIsa::Scalar) return scalar_kernel;
#if FFTLAB_X86_TARGETS
    if (isa == KernelIsa::Avx2) return avx2_kernel;
    if (isa == KernelIsa::Avx512) return avx512_kernel;
#endif
    throw std::invalid_argument("requested FFT kernel is not compiled for this target");
}

KernelIsa tune_kernel(std::size_t n, const std::vector<Swap>& swaps,
                      const std::vector<std::size_t>& offsets, const Vector& twiddles,
                      double& elapsed_ns) {
    const auto capabilities = kernel_capabilities();
    const std::array<KernelIsa, 3> candidates{KernelIsa::Scalar, KernelIsa::Avx2, KernelIsa::Avx512};
    std::array<std::vector<double>, 3> samples;

    Vector seed(n);
    for (std::size_t i = 0; i < n; ++i) {
        seed[i] = {static_cast<double>((17 * i + 3) % 101) / 101.0 - 0.5,
                   static_cast<double>((29 * i + 7) % 103) / 103.0 - 0.5};
    }

    const std::size_t iterations = std::clamp<std::size_t>(262144 / n, 1, 4096);
    const auto tune_start = Clock::now();

    for (std::size_t round = 0; round < 5; ++round) {
        std::array<std::size_t, 3> order{0, 1, 2};
        std::rotate(order.begin(), order.begin() + static_cast<std::ptrdiff_t>(round % 3), order.end());
        if (round & 1U) std::reverse(order.begin(), order.end());

        for (const auto index : order) {
            const auto isa = candidates[index];
            if ((isa == KernelIsa::Avx2 && !capabilities.avx2) ||
                (isa == KernelIsa::Avx512 && !capabilities.avx512)) {
                continue;
            }

            auto data = seed;
            const auto fn = function_for(isa);
            const auto start = Clock::now();
            for (std::size_t iteration = 0; iteration < iterations; ++iteration) {
                fn(data.data(), n, swaps.data(), swaps.size(), offsets.data(), offsets.size(),
                   twiddles.data(), false);
                fn(data.data(), n, swaps.data(), swaps.size(), offsets.data(), offsets.size(),
                   twiddles.data(), true);
            }
            samples[index].push_back(
                std::chrono::duration<double, std::nano>(Clock::now() - start).count() /
                static_cast<double>(2 * iterations));
        }
    }

    KernelIsa best = KernelIsa::Scalar;
    double best_median = std::numeric_limits<double>::infinity();
    for (std::size_t index = 0; index < candidates.size(); ++index) {
        if (samples[index].empty()) continue;
        auto values = samples[index];
        std::sort(values.begin(), values.end());
        const double median = values[values.size() / 2];
        if (median < best_median) {
            best_median = median;
            best = candidates[index];
        }
    }

    elapsed_ns = std::chrono::duration<double, std::nano>(Clock::now() - tune_start).count();
    return best;
}

} // namespace

std::string_view kernel_name(KernelIsa isa) {
    switch (isa) {
        case KernelIsa::Scalar: return "scalar";
        case KernelIsa::Avx2: return "avx2";
        case KernelIsa::Avx512: return "avx512";
        case KernelIsa::Auto: return "auto";
    }
    return "?";
}

KernelCapabilities kernel_capabilities() noexcept {
    KernelCapabilities out{};
#if FFTLAB_X86_TARGETS
    __builtin_cpu_init();
    out.avx2 = __builtin_cpu_supports("avx2") && __builtin_cpu_supports("fma");
    out.avx512 = __builtin_cpu_supports("avx512f") && __builtin_cpu_supports("avx512dq") &&
                 __builtin_cpu_supports("avx512vl") && __builtin_cpu_supports("fma");
#endif
    return out;
}

KernelRadix2Plan::KernelRadix2Plan(std::size_t n, KernelIsa requested) : n_(n) {
    if (!pow2(n)) throw std::invalid_argument("KernelRadix2Plan requires power-of-two N");

    std::size_t bits = 0;
    for (auto q = n; q > 1; q >>= 1) ++bits;
    for (std::size_t i = 0; i < n; ++i) {
        auto x = i;
        std::size_t reversed = 0;
        for (std::size_t bit = 0; bit < bits; ++bit) {
            reversed = (reversed << 1) | (x & 1U);
            x >>= 1;
        }
        if (i < reversed) swaps_.emplace_back(i, reversed);
    }

    std::size_t total_twiddles = 0;
    for (std::size_t len = 2; len <= n;) {
        stage_offsets_.push_back(total_twiddles);
        total_twiddles += len / 2;
        if (len == n) break;
        len <<= 1;
    }
    twiddles_.resize(total_twiddles);
    for (std::size_t stage = 0, len = 2; stage < stage_offsets_.size(); ++stage, len <<= 1) {
        for (std::size_t j = 0; j < len / 2; ++j) {
            const double angle = -2.0 * pi * static_cast<double>(j) / static_cast<double>(len);
            twiddles_[stage_offsets_[stage] + j] = {std::cos(angle), std::sin(angle)};
        }
    }

    const auto capabilities = kernel_capabilities();
    if (requested == KernelIsa::Avx2 && !capabilities.avx2)
        throw std::invalid_argument("AVX2/FMA kernel unavailable on this CPU/compiler");
    if (requested == KernelIsa::Avx512 && !capabilities.avx512)
        throw std::invalid_argument("AVX-512/FMA kernel unavailable on this CPU/compiler");

    selected_ = requested == KernelIsa::Auto
                    ? tune_kernel(n_, swaps_, stage_offsets_, twiddles_, tuning_ns_)
                    : requested;
    fn_ = function_for(selected_);
}

void KernelRadix2Plan::execute(Vector& data, bool inverse) const {
    if (data.size() != n_) throw std::invalid_argument("KernelRadix2Plan input size mismatch");
    fn_(data.data(), n_, swaps_.data(), swaps_.size(), stage_offsets_.data(),
        stage_offsets_.size(), twiddles_.data(), inverse);
}

void KernelRadix2Plan::forward_inplace(Vector& data) const { execute(data, false); }
void KernelRadix2Plan::inverse_inplace(Vector& data) const { execute(data, true); }

void kernel_tests() {
    std::size_t checks = 0;
    auto require = [&](bool condition, const char* message) {
        ++checks;
        if (!condition) throw std::runtime_error(message);
    };

    const auto capabilities = kernel_capabilities();
    for (std::size_t n : {2U, 4U, 8U, 16U, 64U, 256U, 1024U}) {
        Vector input(n);
        for (std::size_t i = 0; i < n; ++i) {
            input[i] = {static_cast<double>((13 * i) % 37) / 37.0 - 0.5,
                        static_cast<double>((7 * i) % 41) / 41.0 - 0.5};
        }

        auto reference = input;
        radix2_inplace(reference);

        std::vector<KernelIsa> modes{KernelIsa::Scalar, KernelIsa::Auto};
        if (capabilities.avx2) modes.push_back(KernelIsa::Avx2);
        if (capabilities.avx512) modes.push_back(KernelIsa::Avx512);

        for (const auto isa : modes) {
            KernelRadix2Plan plan(n, isa);
            auto actual = input;
            plan.forward_inplace(actual);
            for (std::size_t i = 0; i < n; ++i) {
                require(std::abs(actual[i] - reference[i]) < 2e-12 * (1.0 + std::abs(reference[i])),
                        "kernel forward mismatch");
            }
            plan.inverse_inplace(actual);
            for (std::size_t i = 0; i < n; ++i) {
                require(std::abs(actual[i] - input[i]) < 2e-12 * (1.0 + std::abs(input[i])),
                        "kernel round-trip mismatch");
            }
            require(plan.selected_isa() != KernelIsa::Auto,
                    "auto kernel must resolve to a concrete ISA");
        }
    }

    bool threw = false;
    try {
        KernelRadix2Plan invalid(12);
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    require(threw, "kernel plan must reject non-power-of-two sizes");

    std::cout << "PASS: " << checks << " kernel checks\n";
}

} // namespace fftlab
