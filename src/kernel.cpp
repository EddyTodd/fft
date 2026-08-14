#include "fftlab/kernel.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <stdexcept>
#include <utility>

#if (defined(__GNUC__) || defined(__clang__)) && (defined(__x86_64__) || defined(__i386__))
#include <immintrin.h>
#define FFTLAB_X86_TARGETS 1
#else
#define FFTLAB_X86_TARGETS 0
#endif

namespace fftlab {
namespace {

using Swap = std::pair<std::size_t, std::size_t>;
using KernelFn = void (*)(Complex64*, std::size_t, const Swap*, std::size_t,
                          const std::size_t*, std::size_t, const Complex64*, bool);

void apply_permutation(Complex64* data, const Swap* swaps, std::size_t count) {
    for (std::size_t i = 0; i < count; ++i) {
        std::swap(data[swaps[i].first], data[swaps[i].second]);
    }
}

void scalar_kernel(Complex64* data, std::size_t n, const Swap* swaps, std::size_t swap_count,
                   const std::size_t* offsets, std::size_t stage_count,
                   const Complex64* twiddles, bool inverse) {
    apply_permutation(data, swaps, swap_count);
    std::size_t stage = 0;
    for (std::size_t length = 2; stage < stage_count; ++stage, length <<= 1) {
        const auto* stage_twiddles = twiddles + offsets[stage];
        const auto half = length / 2;
        for (std::size_t base = 0; base < n; base += length) {
            for (std::size_t offset = 0; offset < half; ++offset) {
                const Complex64 twiddle =
                    inverse ? std::conj(stage_twiddles[offset]) : stage_twiddles[offset];
                const auto upper = data[base + offset];
                const auto lower = data[base + offset + half] * twiddle;
                data[base + offset] = upper + lower;
                data[base + offset + half] = upper - lower;
            }
        }
    }

    if (inverse) {
        const double scale = 1.0 / static_cast<double>(n);
        for (std::size_t i = 0; i < n; ++i) {
            data[i] *= scale;
        }
    }
}

#if FFTLAB_X86_TARGETS
__attribute__((target("avx2,fma")))
void avx2_kernel(Complex64* data, std::size_t n, const Swap* swaps, std::size_t swap_count,
                 const std::size_t* offsets, std::size_t stage_count,
                 const Complex64* twiddles, bool inverse) {
    apply_permutation(data, swaps, swap_count);
    const __m256d conjugate_mask = _mm256_setr_pd(1.0, -1.0, 1.0, -1.0);
    std::size_t stage = 0;
    for (std::size_t length = 2; stage < stage_count; ++stage, length <<= 1) {
        const auto* stage_twiddles = twiddles + offsets[stage];
        const auto half = length / 2;
        for (std::size_t base = 0; base < n; base += length) {
            std::size_t offset = 0;
            for (; offset + 2 <= half; offset += 2) {
                const auto upper =
                    _mm256_loadu_pd(reinterpret_cast<const double*>(data + base + offset));
                const auto lower = _mm256_loadu_pd(
                    reinterpret_cast<const double*>(data + base + offset + half));
                auto twiddle = _mm256_loadu_pd(
                    reinterpret_cast<const double*>(stage_twiddles + offset));
                if (inverse) {
                    twiddle = _mm256_mul_pd(twiddle, conjugate_mask);
                }

                const auto twiddle_real = _mm256_movedup_pd(twiddle);
                const auto twiddle_imag = _mm256_permute_pd(twiddle, 0xF);
                const auto swapped_lower = _mm256_permute_pd(lower, 0x5);
                const auto cross = _mm256_mul_pd(swapped_lower, twiddle_imag);
                const auto product = _mm256_fmaddsub_pd(lower, twiddle_real, cross);
                _mm256_storeu_pd(reinterpret_cast<double*>(data + base + offset),
                                 _mm256_add_pd(upper, product));
                _mm256_storeu_pd(reinterpret_cast<double*>(data + base + offset + half),
                                 _mm256_sub_pd(upper, product));
            }

            for (; offset < half; ++offset) {
                const Complex64 twiddle =
                    inverse ? std::conj(stage_twiddles[offset]) : stage_twiddles[offset];
                const auto upper = data[base + offset];
                const auto lower = data[base + offset + half] * twiddle;
                data[base + offset] = upper + lower;
                data[base + offset + half] = upper - lower;
            }
        }
    }

    if (inverse) {
        const auto scale = _mm256_set1_pd(1.0 / static_cast<double>(n));
        std::size_t i = 0;
        for (; i + 2 <= n; i += 2) {
            const auto value = _mm256_loadu_pd(reinterpret_cast<const double*>(data + i));
            _mm256_storeu_pd(reinterpret_cast<double*>(data + i), _mm256_mul_pd(value, scale));
        }
        for (; i < n; ++i) {
            data[i] /= static_cast<double>(n);
        }
    }
}

__attribute__((target("avx512f,avx512dq,avx512vl,fma")))
void avx512_kernel(Complex64* data, std::size_t n, const Swap* swaps, std::size_t swap_count,
                   const std::size_t* offsets, std::size_t stage_count,
                   const Complex64* twiddles, bool inverse) {
    apply_permutation(data, swaps, swap_count);
    const __m512d multiply_sign =
        _mm512_setr_pd(-1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0);
    const __m512d conjugate_mask =
        _mm512_setr_pd(1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0);
    std::size_t stage = 0;
    for (std::size_t length = 2; stage < stage_count; ++stage, length <<= 1) {
        const auto* stage_twiddles = twiddles + offsets[stage];
        const auto half = length / 2;
        for (std::size_t base = 0; base < n; base += length) {
            std::size_t offset = 0;
            for (; offset + 4 <= half; offset += 4) {
                const auto upper =
                    _mm512_loadu_pd(reinterpret_cast<const void*>(data + base + offset));
                const auto lower = _mm512_loadu_pd(
                    reinterpret_cast<const void*>(data + base + offset + half));
                auto twiddle =
                    _mm512_loadu_pd(reinterpret_cast<const void*>(stage_twiddles + offset));
                if (inverse) {
                    twiddle = _mm512_mul_pd(twiddle, conjugate_mask);
                }

                const auto twiddle_real = _mm512_movedup_pd(twiddle);
                const auto twiddle_imag = _mm512_permute_pd(twiddle, 0xFF);
                const auto swapped_lower = _mm512_permute_pd(lower, 0x55);
                const auto cross =
                    _mm512_mul_pd(_mm512_mul_pd(swapped_lower, twiddle_imag), multiply_sign);
                const auto product = _mm512_fmadd_pd(lower, twiddle_real, cross);
                _mm512_storeu_pd(reinterpret_cast<void*>(data + base + offset),
                                 _mm512_add_pd(upper, product));
                _mm512_storeu_pd(reinterpret_cast<void*>(data + base + offset + half),
                                 _mm512_sub_pd(upper, product));
            }

            for (; offset < half; ++offset) {
                const Complex64 twiddle =
                    inverse ? std::conj(stage_twiddles[offset]) : stage_twiddles[offset];
                const auto upper = data[base + offset];
                const auto lower = data[base + offset + half] * twiddle;
                data[base + offset] = upper + lower;
                data[base + offset + half] = upper - lower;
            }
        }
    }

    if (inverse) {
        const auto scale = _mm512_set1_pd(1.0 / static_cast<double>(n));
        std::size_t i = 0;
        for (; i + 4 <= n; i += 4) {
            const auto value = _mm512_loadu_pd(reinterpret_cast<const void*>(data + i));
            _mm512_storeu_pd(reinterpret_cast<void*>(data + i), _mm512_mul_pd(value, scale));
        }
        for (; i < n; ++i) {
            data[i] /= static_cast<double>(n);
        }
    }
}
#endif

[[nodiscard]] KernelFn function_for(KernelIsa isa) {
    if (isa == KernelIsa::Scalar) {
        return scalar_kernel;
    }
#if FFTLAB_X86_TARGETS
    if (isa == KernelIsa::Avx2) {
        return avx2_kernel;
    }
    if (isa == KernelIsa::Avx512) {
        return avx512_kernel;
    }
#endif
    throw std::invalid_argument("requested FFT kernel is not compiled for this target");
}

[[nodiscard]] KernelCapabilities detect_kernel_capabilities() noexcept {
    KernelCapabilities capabilities{};
#if FFTLAB_X86_TARGETS
    __builtin_cpu_init();
    capabilities.avx2 =
        __builtin_cpu_supports("avx2") && __builtin_cpu_supports("fma");
    capabilities.avx512 = __builtin_cpu_supports("avx512f") &&
                          __builtin_cpu_supports("avx512dq") &&
                          __builtin_cpu_supports("avx512vl") &&
                          __builtin_cpu_supports("fma");
#endif
    return capabilities;
}

}  // namespace

std::string_view kernel_name(KernelIsa isa) noexcept {
    switch (isa) {
        case KernelIsa::Scalar:
            return "scalar";
        case KernelIsa::Avx2:
            return "avx2";
        case KernelIsa::Avx512:
            return "avx512";
    }
    return "unknown";
}

KernelCapabilities kernel_capabilities() noexcept {
    static const KernelCapabilities capabilities = detect_kernel_capabilities();
    return capabilities;
}

KernelRadix2Plan::KernelRadix2Plan(std::size_t n, KernelIsa requested)
    : n_(n), selected_(requested) {
    if (n == 0 || !pow2(n)) {
        throw std::invalid_argument("KernelRadix2Plan requires power-of-two N >= 1");
    }

    const auto bits = ilog2(n);
    swaps_.reserve(n / 2);
    for (std::size_t i = 0; i < n; ++i) {
        auto value = i;
        std::size_t reversed = 0;
        for (std::size_t bit = 0; bit < bits; ++bit) {
            reversed = (reversed << 1) | (value & 1U);
            value >>= 1;
        }
        if (i < reversed) {
            swaps_.emplace_back(i, reversed);
        }
    }

    stage_offsets_.reserve(bits);
    std::size_t total_twiddles = 0;
    for (std::size_t length = 2; length <= n;) {
        stage_offsets_.push_back(total_twiddles);
        total_twiddles += length / 2;
        if (length == n) {
            break;
        }
        length <<= 1;
    }

    twiddles_.resize(total_twiddles);
    for (std::size_t stage = 0, length = 2; stage < stage_offsets_.size();
         ++stage, length <<= 1) {
        for (std::size_t offset = 0; offset < length / 2; ++offset) {
            const double angle = -2.0 * std::numbers::pi_v<double> *
                                 static_cast<double>(offset) / static_cast<double>(length);
            twiddles_[stage_offsets_[stage] + offset] = {std::cos(angle), std::sin(angle)};
        }
    }

    const auto capabilities = kernel_capabilities();
    if (requested == KernelIsa::Avx2 && !capabilities.avx2) {
        throw std::invalid_argument("AVX2/FMA kernel unavailable");
    }
    if (requested == KernelIsa::Avx512 && !capabilities.avx512) {
        throw std::invalid_argument("AVX-512/FMA kernel unavailable");
    }
    fn_ = function_for(requested);
}

void KernelRadix2Plan::execute(std::span<Complex64> data, bool inverse) const {
    if (data.size() != n_) {
        throw std::invalid_argument("KernelRadix2Plan buffer size mismatch");
    }
    if (n_ <= 1) {
        return;
    }
    fn_(data.data(), n_, swaps_.data(), swaps_.size(), stage_offsets_.data(),
        stage_offsets_.size(), twiddles_.data(), inverse);
}

void KernelRadix2Plan::forward_inplace(std::span<Complex64> data) const {
    execute(data, false);
}

void KernelRadix2Plan::inverse_inplace(std::span<Complex64> data) const {
    execute(data, true);
}

}  // namespace fftlab
