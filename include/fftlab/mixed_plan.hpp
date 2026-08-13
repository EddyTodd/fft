#pragma once

#include "fftlab/codelet.hpp"

namespace fftlab {

template <FftScalar T = double>
class MixedRadixPlan {
    struct Node {
        std::size_t n{};
        std::size_t radix{};
        std::size_t m{};
        bool leaf{};
        VectorT<T> twiddles;
        std::unique_ptr<SmallDftCodelet<T>> codelet;
        std::unique_ptr<Node> child;
    };

public:
    explicit MixedRadixPlan(std::size_t n) : n_(n), root_(build(n)) {
        if (n == 0) throw std::invalid_argument("MixedRadixPlan requires N >= 1");
    }
    MixedRadixPlan(MixedRadixPlan&&) noexcept = default;
    MixedRadixPlan& operator=(MixedRadixPlan&&) noexcept = default;
    MixedRadixPlan(const MixedRadixPlan&) = delete;
    MixedRadixPlan& operator=(const MixedRadixPlan&) = delete;

    [[nodiscard]] std::size_t size() const noexcept { return n_; }
    [[nodiscard]] std::size_t scratch_size() const noexcept { return n_; }
    [[nodiscard]] std::size_t top_radix() const noexcept { return root_->radix; }

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
    [[nodiscard]] static std::unique_ptr<Node> build(std::size_t n) {
        if (n == 0) throw std::invalid_argument("MixedRadixPlan requires N >= 1");
        auto node = std::make_unique<Node>();
        node->n = n;
        if (n == 1) { node->leaf = true; node->radix = 1; node->m = 1; return node; }
        if (small_codelet_supported(n)) { node->leaf = true; node->radix = n; node->m = 1; node->codelet = std::make_unique<SmallDftCodelet<T>>(n); return node; }
        const auto radix = choose_small_radix(n);
        if (radix == 0) {
            // A non-codelet prime/rough factor remains a planned direct leaf:
            // precompute its complete DFT matrix once so execution performs no
            // allocation or trigonometric setup. Structural top-level planning
            // still prefers Bluestein for rough composites/primes.
            node->leaf = true;
            node->radix = n;
            node->m = 1;
            node->twiddles.resize(n * n);
            for (std::size_t k = 0; k < n; ++k)
                for (std::size_t q = 0; q < n; ++q)
                    node->twiddles[k * n + q] = root<T>(-T{2} * std::numbers::pi_v<T> *
                        static_cast<T>(k * q) / static_cast<T>(n));
            return node;
        }
        node->radix = radix;
        node->m = n / radix;
        node->codelet = std::make_unique<SmallDftCodelet<T>>(radix);
        node->child = build(node->m);
        node->twiddles.resize(node->m * radix);
        for (std::size_t k0 = 0; k0 < node->m; ++k0)
            for (std::size_t q = 0; q < radix; ++q)
                node->twiddles[k0 * radix + q] = root<T>(-T{2} * std::numbers::pi_v<T> *
                    static_cast<T>(k0 * q) / static_cast<T>(n));
        return node;
    }

    static void direct_leaf(const Node& node, std::span<ComplexT<T>> data,
                            std::span<ComplexT<T>> scratch, Direction direction) {
        const auto n = node.n;
        for (std::size_t k = 0; k < n; ++k) {
            ComplexT<T> sum{};
            for (std::size_t q = 0; q < n; ++q) {
                auto w = node.twiddles[k * n + q];
                if (direction == Direction::Inverse) w = std::conj(w);
                sum += data[q] * w;
            }
            scratch[k] = sum;
        }
        if (direction == Direction::Inverse) {
            const auto scale = T{1} / static_cast<T>(n);
            for (std::size_t k = 0; k < n; ++k) scratch[k] *= scale;
        }
        std::copy(scratch.begin(), scratch.begin() + static_cast<std::ptrdiff_t>(n), data.begin());
    }

    static void execute_node(const Node& node, std::span<ComplexT<T>> data,
                             std::span<ComplexT<T>> scratch, Direction direction) {
        if (node.n == 1) return;
        if (node.leaf) {
            if (node.codelet) node.codelet->execute(data.first(node.n), direction);
            else direct_leaf(node, data.first(node.n), scratch.first(node.n), direction);
            return;
        }
        const auto r = node.radix, m = node.m;
        for (std::size_t q = 0; q < r; ++q)
            for (std::size_t j = 0; j < m; ++j)
                scratch[q * m + j] = data[r * j + q];
        for (std::size_t q = 0; q < r; ++q)
            execute_node(*node.child, scratch.subspan(q * m, m), data.subspan(q * m, m), direction);
        std::array<ComplexT<T>, 7> values{};
        for (std::size_t k0 = 0; k0 < m; ++k0) {
            for (std::size_t q = 0; q < r; ++q) {
                auto w = node.twiddles[k0 * r + q];
                if (direction == Direction::Inverse) w = std::conj(w);
                values[q] = scratch[q * m + k0] * w;
            }
            node.codelet->execute(std::span<ComplexT<T>>(values.data(), r), direction);
            for (std::size_t k1 = 0; k1 < r; ++k1) data[k0 + m * k1] = values[k1];
        }
    }

    void execute(std::span<ComplexT<T>> data, std::span<ComplexT<T>> scratch, Direction direction) const {
        if (data.size() != n_ || scratch.size() < n_) throw std::invalid_argument("MixedRadixPlan buffer size mismatch");
        execute_node(*root_, data, scratch.first(n_), direction);
    }
    void copy_checked(std::span<const ComplexT<T>> input, std::span<ComplexT<T>> output,
                      std::span<ComplexT<T>> scratch) const {
        if (input.size() != n_ || output.size() != n_ || scratch.size() < n_)
            throw std::invalid_argument("MixedRadixPlan buffer size mismatch");
        std::copy(input.begin(), input.end(), output.begin());
    }

    std::size_t n_{};
    std::unique_ptr<Node> root_;
};
MixedRadixPlan(std::size_t) -> MixedRadixPlan<double>;

} // namespace fftlab
