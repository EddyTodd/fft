#include "fftlab/fft.hpp"

#include <charconv>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace fftlab {
struct Opt {
    bool help = false, self = false, verify = false, accuracy_suite = false, complexity = false, one = false, suite = false, csv = false, raw_csv = false;
    Algo algo = Algo::Auto; SignalKind sig = SignalKind::Tones; std::size_t size = 1024, samples = 31, warm = 5; double target = 5;
    std::vector<std::size_t> sizes{8, 16, 32, 64, 128, 256, 512, 1000, 1009, 1024, 2048, 4096};
};
std::size_t number(std::string_view s) { std::size_t n{}; auto [p, e] = std::from_chars(s.data(), s.data() + s.size(), n); if (e != std::errc{} || p != s.data() + s.size() || !n) throw std::invalid_argument("invalid integer"); return n; }
std::vector<std::size_t> list(std::string_view s) { std::vector<std::size_t> v; while (!s.empty()) { const auto p = s.find(','); v.push_back(number(s.substr(0, p))); if (p == s.npos) break; s.remove_prefix(p + 1); } return v; }
Opt options(int ac, char** av) {
    Opt o; auto val = [&](int& i) { if (++i >= ac) throw std::invalid_argument("missing option value"); return std::string_view(av[i]); };
    for (int i = 1; i < ac; ++i) {
        const std::string_view a = av[i];
        if (a == "-h" || a == "--help") o.help = true; else if (a == "--self-test") o.self = true; else if (a == "--verify") o.verify = true;
        else if (a == "--accuracy-suite") o.accuracy_suite = true; else if (a == "--complexity") o.complexity = true; else if (a == "--benchmark") o.one = true;
        else if (a == "--benchmark-suite") o.suite = true; else if (a == "--csv") o.csv = true; else if (a == "--raw-csv") o.raw_csv = true;
        else if (a == "--algorithm") o.algo = parse_algo(val(i)); else if (a == "--signal") o.sig = parse_signal(val(i)); else if (a == "--size") o.size = number(val(i));
        else if (a == "--samples") o.samples = number(val(i)); else if (a == "--warmups") o.warm = number(val(i)); else if (a == "--target-ms") o.target = std::stod(std::string(val(i)));
        else if (a == "--sizes") o.sizes = list(val(i)); else throw std::invalid_argument("unknown option");
    }
    return o;
}
void help() {
    std::cout << "fft - dependency-free C++23 Fourier transform research laboratory\n\n"
              << "--self-test\n--verify --algorithm ALGO --size N [--signal KIND]\n--accuracy-suite [--sizes ...] [--csv]\n"
              << "--complexity --algorithm ALGO --size N [--csv]\n--benchmark --algorithm ALGO --size N [--samples 31] [--target-ms 5] [--csv|--raw-csv]\n"
              << "--benchmark-suite [--sizes ...] [--csv|--raw-csv]\n\nAlgorithms: auto, dft, radix2-iterative, radix2-recursive, stockham-radix2, radix4, split-radix, mixed-radix, rader, bluestein\n"
              << "Signals: random, tones, impulse, alternating, dynamic-range\n";
}
void print_stats(const Stats& s, bool csv) {
    if (csv) std::cout << name(s.algo) << ',' << s.n << ',' << s.samples << ',' << s.iters << ',' << std::setprecision(12) << s.min << ',' << s.p05 << ',' << s.median << ',' << s.mean << ',' << s.p95 << ',' << s.max << ',' << s.sd << ',' << s.mad << ',' << s.ci_lo << ',' << s.ci_hi << ',' << s.equiv_gflops << '\n';
    else std::cout << std::left << std::setw(20) << name(s.algo) << " N=" << std::setw(6) << s.n << " median=" << std::fixed << std::setprecision(1) << s.median << " ns  95% bootstrap CI=[" << s.ci_lo << ", " << s.ci_hi << "]  MAD=" << s.mad << "  radix2-equiv=" << std::setprecision(3) << s.equiv_gflops << " GFLOP/s\n";
}
void print_raw(const Stats& s) { for (std::size_t i = 0; i < s.raw.size(); ++i) std::cout << name(s.algo) << ',' << s.n << ',' << i << ',' << s.iters << ',' << std::setprecision(12) << s.raw[i] << '\n'; }
void print_accuracy(const Accuracy& a, bool csv) {
    if (csv) std::cout << name(a.algo) << ',' << a.n << ',' << signal_name(a.sig) << ',' << std::setprecision(12) << a.forward.l1 << ',' << a.forward.l2 << ',' << a.forward.linf << ',' << a.backward.l1 << ',' << a.backward.l2 << ',' << a.backward.linf << ',' << a.roundtrip_max << '\n';
    else std::cout << name(a.algo) << " N=" << a.n << " signal=" << signal_name(a.sig) << " forward(L1,L2,Linf)=(" << a.forward.l1 << ',' << a.forward.l2 << ',' << a.forward.linf << ") backward=(" << a.backward.l1 << ',' << a.backward.l2 << ',' << a.backward.linf << ") roundtrip_max=" << a.roundtrip_max << '\n';
}

} // namespace fftlab

int main(int argc, char** argv) {
    using namespace fftlab;
    try {
        const auto o = options(argc, argv); if (argc == 1 || o.help) { help(); return 0; } if (o.self) { tests(); return 0; }
        if (o.verify) {
            if (o.size > 4096) throw std::invalid_argument("verification capped at N=4096");
            print_accuracy(accuracy(signal(o.size, o.sig), o.algo, o.sig), false); return 0;
        }
        if (o.accuracy_suite) {
            if (o.csv) std::cout << "algorithm,N,signal,forward_l1,forward_l2,forward_linf,backward_l1,backward_l2,backward_linf,roundtrip_max_abs\n";
            for (auto n : o.sizes) {
                if (n > 2048) continue;
                for (auto sk : signal_kinds) { const auto x = signal(n, sk); for (auto a : suite(n)) print_accuracy(accuracy(x, a, sk), o.csv); }
            }
            return 0;
        }
        if (o.complexity) {
            const auto m = model(o.algo, o.size);
            if (o.csv) std::cout << "algorithm,N,structural_complex_adds,structural_complex_multiplies,peak_workspace_complex,note\n" << name(o.algo) << ',' << o.size << ',' << m.adds << ',' << m.muls << ',' << m.workspace << ',' << '"' << m.note << '"' << '\n';
            else std::cout << "algorithm: " << name(o.algo) << "\nN: " << o.size << "\nstructural_complex_adds: " << m.adds << "\nstructural_complex_multiplies: " << m.muls << "\npeak_workspace_complex_estimate: " << m.workspace << "\nmodel_note: " << m.note << '\n';
            return 0;
        }
        if (o.one) {
            const auto s = bench(signal(o.size, o.sig), o.algo, o.samples, o.warm, o.target);
            if (o.raw_csv) { std::cout << "algorithm,N,sample,iterations_per_sample,ns_per_transform\n"; print_raw(s); }
            else { if (o.csv) std::cout << "algorithm,N,samples,iterations_per_sample,min_ns,p05_ns,median_ns,mean_ns,p95_ns,max_ns,stddev_ns,mad_ns,median_ci95_low_ns,median_ci95_high_ns,radix2_equiv_gflops\n"; print_stats(s, o.csv); }
            return 0;
        }
        if (o.suite) {
            if (o.raw_csv) std::cout << "algorithm,N,sample,iterations_per_sample,ns_per_transform\n";
            else if (o.csv) std::cout << "algorithm,N,samples,iterations_per_sample,min_ns,p05_ns,median_ns,mean_ns,p95_ns,max_ns,stddev_ns,mad_ns,median_ci95_low_ns,median_ci95_high_ns,radix2_equiv_gflops\n";
            for (auto n : o.sizes) for (auto a : suite(n)) { const auto s = bench(signal(n, o.sig), a, o.samples, o.warm, o.target); if (o.raw_csv) print_raw(s); else print_stats(s, o.csv); }
            return 0;
        }
        throw std::invalid_argument("no action requested");
    } catch (const std::exception& e) { std::cerr << "error: " << e.what() << '\n'; return 2; }
}
