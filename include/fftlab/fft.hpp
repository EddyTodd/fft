#pragma once

#include <array>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

namespace fftlab {
using Complex = std::complex<double>;
using Vector = std::vector<Complex>;

enum class Algo { Auto, Dft, Radix2, Recursive, Stockham, Radix4, SplitRadix, Mixed, Rader, Bluestein };
inline constexpr std::array<Algo, 10> all_algos{Algo::Auto, Algo::Dft, Algo::Radix2, Algo::Recursive, Algo::Stockham, Algo::Radix4, Algo::SplitRadix, Algo::Mixed, Algo::Rader, Algo::Bluestein};

enum class SignalKind { Random, Tones, Impulse, Alternating, DynamicRange };
inline constexpr std::array<SignalKind, 5> signal_kinds{SignalKind::Random, SignalKind::Tones, SignalKind::Impulse, SignalKind::Alternating, SignalKind::DynamicRange};

struct NormErr { double l1{}, l2{}, linf{}; };
struct Accuracy { Algo algo{}; std::size_t n{}; SignalKind sig{}; NormErr forward{}, backward{}; double roundtrip_max{}; };
struct Model { long double adds{}, muls{}, workspace{}; std::string note; };
struct Stats {
    std::size_t n{}, samples{}, iters{}; Algo algo{};
    double min{}, p05{}, median{}, mean{}, p95{}, max{}, sd{}, mad{}, ci_lo{}, ci_hi{}, equiv_gflops{};
    std::vector<double> raw;
};

bool pow2(std::size_t n);
bool is_prime(std::size_t n);
std::string_view name(Algo a);
Algo parse_algo(std::string_view s);
bool supports(Algo a, std::size_t n);
std::string_view signal_name(SignalKind s);
SignalKind parse_signal(std::string_view s);

Vector dft(const Vector& x, bool inv = false);
void radix2_inplace(Vector& a, bool inv = false);
Vector radix2(const Vector& x, bool inv = false);
Vector recursive(const Vector& x, bool inv = false);
Vector stockham(const Vector& x, bool inv = false);
Vector radix4(const Vector& x, bool inv = false);
Vector split_radix(const Vector& x, bool inv = false);
Vector mixed(const Vector& x, bool inv = false);
Vector bluestein(const Vector& x, bool inv = false);
Vector rader(const Vector& x, bool inv = false);
Vector transform(const Vector& x, Algo a, bool inv = false);

Vector signal(std::size_t n, SignalKind kind = SignalKind::Tones, std::uint64_t seed = 0xF17F17ULL);
Accuracy accuracy(const Vector& x, Algo a, SignalKind sig);
Model model(Algo a, std::size_t n);
Stats bench(const Vector& x, Algo a, std::size_t samples, std::size_t warm, double target);
std::vector<Algo> suite(std::size_t n);
void tests();
}
