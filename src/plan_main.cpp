#include "fftlab/plan.hpp"

#include <charconv>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace {
std::size_t number(std::string_view s) { std::size_t n{}; auto [p, e] = std::from_chars(s.data(), s.data() + s.size(), n); if (e != std::errc{} || p != s.data() + s.size() || !n) throw std::invalid_argument("invalid integer"); return n; }
double real_number(std::string_view s) { return std::stod(std::string(s)); }
std::vector<std::size_t> list(std::string_view s) { std::vector<std::size_t> v; while (!s.empty()) { const auto p = s.find(','); v.push_back(number(s.substr(0, p))); if (p == s.npos) break; s.remove_prefix(p + 1); } return v; }
void print_summary(const fftlab::PlanBenchmark& b, bool csv) {
    if (csv) {
        std::cout << b.n << ',' << b.samples << ',' << b.iterations_per_sample << ',' << std::setprecision(12)
                  << b.complex_setup.median << ',' << b.real_setup.median << ',' << b.legacy_complex.median << ',' << b.planned_complex.median << ',' << b.planned_real.median << ','
                  << b.plan_speedup << ',' << b.real_speedup << ',' << b.complex_setup_break_even_transforms << ',' << b.real_setup_break_even_transforms << ','
                  << b.legacy_complex.ci_lo << ',' << b.legacy_complex.ci_hi << ',' << b.planned_complex.ci_lo << ',' << b.planned_complex.ci_hi << ',' << b.planned_real.ci_lo << ',' << b.planned_real.ci_hi << '\n';
    } else {
        std::cout << "N=" << b.n << " complex_setup=" << b.complex_setup.median << " ns real_setup=" << b.real_setup.median << " ns legacy=" << b.legacy_complex.median << " ns planned=" << b.planned_complex.median
                  << " ns real=" << b.planned_real.median << " ns plan_speedup=" << b.plan_speedup << "x real_speedup=" << b.real_speedup
                  << "x complex_break_even=" << b.complex_setup_break_even_transforms << " transforms real_break_even=" << b.real_setup_break_even_transforms << " transforms\n";
    }
}
void print_raw(const fftlab::PlanBenchmark& b) {
    const auto emit = [&](std::string_view mode, const fftlab::PlanDistribution& d) { for (std::size_t i = 0; i < d.raw.size(); ++i) std::cout << b.n << ',' << mode << ',' << i << ',' << b.iterations_per_sample << ',' << std::setprecision(12) << d.raw[i] << '\n'; };
    emit("complex-setup", b.complex_setup); emit("real-setup", b.real_setup); emit("legacy-complex", b.legacy_complex); emit("planned-complex", b.planned_complex); emit("planned-real", b.planned_real);
}
}
int main(int argc, char** argv) {
    try {
        bool self = false, bench = false, sweep = false, csv = false, raw = false; std::size_t n = 1024, samples = 31, warmups = 5; std::uint64_t seed = 0xF17F17ULL; double target = 5.0; std::vector<std::size_t> sizes{256, 1024, 4096, 16384};
        auto value = [&](int& i) { if (++i >= argc) throw std::invalid_argument("missing option value"); return std::string_view(argv[i]); };
        for (int i = 1; i < argc; ++i) { const std::string_view a = argv[i]; if (a == "--self-test") self = true; else if (a == "--benchmark") bench = true; else if (a == "--sweep") sweep = true; else if (a == "--csv") csv = true; else if (a == "--raw-csv") raw = true; else if (a == "--size") n = number(value(i)); else if (a == "--sizes") sizes = list(value(i)); else if (a == "--samples") samples = number(value(i)); else if (a == "--warmups") warmups = number(value(i)); else if (a == "--target-ms") target = real_number(value(i)); else if (a == "--seed") seed = static_cast<std::uint64_t>(number(value(i))); else if (a == "-h" || a == "--help") { std::cout << "fft-plan --self-test | --benchmark --size N | --sweep [--sizes ...] [--samples 31] [--target-ms 5] [--seed N] [--csv|--raw-csv]\n"; return 0; } else throw std::invalid_argument("unknown option"); }
        if (self) { fftlab::planned_tests(); return 0; }
        if (csv && !raw) std::cout << "N,samples,iterations_per_sample,complex_setup_median_ns,real_setup_median_ns,legacy_complex_median_ns,planned_complex_median_ns,planned_real_median_ns,plan_speedup,real_speedup,complex_setup_break_even_transforms,real_setup_break_even_transforms,legacy_ci95_low_ns,legacy_ci95_high_ns,planned_ci95_low_ns,planned_ci95_high_ns,real_ci95_low_ns,real_ci95_high_ns\n";
        if (raw) std::cout << "N,mode,sample,iterations_per_sample,ns\n";
        if (bench) { const auto b = fftlab::benchmark_plans(n, samples, warmups, target, seed); if (raw) print_raw(b); else print_summary(b, csv); return 0; }
        if (sweep) { for (auto size : sizes) { const auto b = fftlab::benchmark_plans(size, samples, warmups, target, seed ^ size); if (raw) print_raw(b); else print_summary(b, csv); } return 0; }
        throw std::invalid_argument("no action requested");
    } catch (const std::exception& e) { std::cerr << "error: " << e.what() << '\n'; return 2; }
}
