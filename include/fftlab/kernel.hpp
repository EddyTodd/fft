#pragma once

#include "fftlab/fft.hpp"

#include <cstddef>
#include <string_view>
#include <utility>
#include <vector>

namespace fftlab {

enum class KernelIsa { Scalar, Avx2, Avx512, Auto };

struct KernelCapabilities {
    bool avx2{};
    bool avx512{};
};

std::string_view kernel_name(KernelIsa isa);
KernelCapabilities kernel_capabilities() noexcept;

class KernelRadix2Plan {
public:
    explicit KernelRadix2Plan(std::size_t n, KernelIsa requested = KernelIsa::Auto);

    [[nodiscard]] std::size_t size() const noexcept { return n_; }
    [[nodiscard]] KernelIsa selected_isa() const noexcept { return selected_; }
    [[nodiscard]] double tuning_ns() const noexcept { return tuning_ns_; }
    [[nodiscard]] std::size_t stored_twiddles() const noexcept { return twiddles_.size(); }
    [[nodiscard]] std::size_t stored_swaps() const noexcept { return swaps_.size(); }

    void forward_inplace(Vector& data) const;
    void inverse_inplace(Vector& data) const;

private:
    using Swap = std::pair<std::size_t, std::size_t>;
    using KernelFn = void (*)(Complex*, std::size_t, const Swap*, std::size_t,
                              const std::size_t*, std::size_t, const Complex*, bool);

    void execute(Vector& data, bool inverse) const;

    KernelFn fn_{};
    std::size_t n_{};
    KernelIsa selected_{KernelIsa::Scalar};
    double tuning_ns_{};
    std::vector<Swap> swaps_;
    std::vector<std::size_t> stage_offsets_;
    Vector twiddles_;
};

void kernel_tests();

} // namespace fftlab
