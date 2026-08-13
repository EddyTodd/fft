#pragma once

#include "fftlab/fft.hpp"

#include <cstddef>
#include <span>
#include <string_view>
#include <utility>
#include <vector>

namespace fftlab {

// v1 SIMD kernels are intentionally an explicit binary64 extension. The generic
// Radix2Plan<float/double> is the portable planner path; architecture selection
// here is never timed or auto-tuned by the library.
enum class KernelIsa { Scalar, Avx2, Avx512 };

struct KernelCapabilities {
    bool avx2{};
    bool avx512{};
};

[[nodiscard]] std::string_view kernel_name(KernelIsa isa) noexcept;
[[nodiscard]] KernelCapabilities kernel_capabilities() noexcept;

class KernelRadix2Plan {
public:
    explicit KernelRadix2Plan(std::size_t n, KernelIsa requested = KernelIsa::Scalar);

    [[nodiscard]] std::size_t size() const noexcept { return n_; }
    [[nodiscard]] KernelIsa selected_isa() const noexcept { return selected_; }
    [[nodiscard]] std::size_t stored_twiddles() const noexcept { return twiddles_.size(); }
    [[nodiscard]] std::size_t stored_swaps() const noexcept { return swaps_.size(); }
    [[nodiscard]] std::size_t scratch_size() const noexcept { return 0; }

    void forward_inplace(std::span<Complex64> data) const;
    void inverse_inplace(std::span<Complex64> data) const;

private:
    using Swap = std::pair<std::size_t, std::size_t>;
    using KernelFn = void (*)(Complex64*, std::size_t, const Swap*, std::size_t,
                              const std::size_t*, std::size_t, const Complex64*, bool);

    void execute(std::span<Complex64> data, bool inverse) const;

    KernelFn fn_{};
    std::size_t n_{};
    KernelIsa selected_{KernelIsa::Scalar};
    std::vector<Swap> swaps_;
    std::vector<std::size_t> stage_offsets_;
    Vector64 twiddles_;
};

} // namespace fftlab
