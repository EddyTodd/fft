#include "fftlab/arbitrary_plan.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cstddef>
#include <cstdint>
#ifdef _WIN32
#include <windows.h>
#else
#include <dlfcn.h>
#endif
#include <iomanip>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>

namespace {
using Clock = std::chrono::steady_clock;
using FftwPlan = void*;
constexpr int fftw_forward = -1;
constexpr int fftw_backward = 1;
constexpr unsigned fftw_measure = 0U;
constexpr unsigned fftw_estimate = 1U << 6;

#ifdef _WIN32
using LibraryHandle = HMODULE;
LibraryHandle open_library(const char* name) { return LoadLibraryA(name); }
void* load_symbol(LibraryHandle handle, const char* name) { return reinterpret_cast<void*>(GetProcAddress(handle, name)); }
void close_library(LibraryHandle handle) { if (handle) FreeLibrary(handle); }
#else
using LibraryHandle = void*;
LibraryHandle open_library(const char* name) { return dlopen(name, RTLD_NOW | RTLD_LOCAL); }
void* load_symbol(LibraryHandle handle, const char* name) { return dlsym(handle, name); }
void close_library(LibraryHandle handle) { if (handle) dlclose(handle); }
#endif

struct FftwApi {
    using PlanDft = FftwPlan (*)(int, double*, double*, int, unsigned);
    using ExecuteDft = void (*)(FftwPlan, double*, double*);
    using Destroy = void (*)(FftwPlan);
    using Malloc = void* (*)(std::size_t);
    using Free = void (*)(void*);
    using ForgetWisdom = void (*)();
    using Cleanup = void (*)();

    LibraryHandle handle{};
    std::string library;
    std::string version;
    PlanDft plan_dft{};
    ExecuteDft execute_dft{};
    Destroy destroy_plan{};
    Malloc malloc_fn{};
    Free free_fn{};
    ForgetWisdom forget_wisdom{};
    Cleanup cleanup{};

    FftwApi() {
        for (const char* candidate : {"libfftw3.so.3", "libfftw3.so", "libfftw3.3.dylib",
                                      "libfftw3.dylib", "libfftw3-3.dll", "fftw3.dll"}) {
            handle = open_library(candidate);
            if (handle) { library = candidate; break; }
        }
        if (!handle) return;
        auto symbol = [&](auto& function, const char* name) {
            function = reinterpret_cast<std::remove_reference_t<decltype(function)>>(load_symbol(handle, name));
            if (!function) throw std::runtime_error(std::string("missing FFTW symbol: ") + name);
        };
        symbol(plan_dft, "fftw_plan_dft_1d");
        symbol(execute_dft, "fftw_execute_dft");
        symbol(destroy_plan, "fftw_destroy_plan");
        symbol(malloc_fn, "fftw_malloc");
        symbol(free_fn, "fftw_free");
        forget_wisdom = reinterpret_cast<ForgetWisdom>(load_symbol(handle, "fftw_forget_wisdom"));
        cleanup = reinterpret_cast<Cleanup>(load_symbol(handle, "fftw_cleanup"));
        if (const auto* text = static_cast<const char*>(load_symbol(handle, "fftw_version"))) version = text;
    }
    ~FftwApi() { if (cleanup) cleanup(); close_library(handle); }
    FftwApi(const FftwApi&) = delete;
    FftwApi& operator=(const FftwApi&) = delete;
    explicit operator bool() const noexcept { return handle != nullptr; }
};

struct FftwBuffer {
    FftwApi* api{};
    double* data{};
    FftwBuffer(FftwApi& owner, std::size_t doubles)
        : api(&owner), data(static_cast<double*>(owner.malloc_fn(sizeof(double) * doubles))) {
        if (!data) throw std::bad_alloc();
    }
    ~FftwBuffer() { api->free_fn(data); }
    FftwBuffer(const FftwBuffer&) = delete;
    FftwBuffer& operator=(const FftwBuffer&) = delete;
};

struct FftwPlanHandle {
    FftwApi* api{};
    FftwPlan plan{};
    FftwPlanHandle(FftwApi& owner, FftwPlan value) : api(&owner), plan(value) {
        if (!plan) throw std::runtime_error("FFTW plan creation failed");
    }
    ~FftwPlanHandle() { api->destroy_plan(plan); }
    FftwPlanHandle(const FftwPlanHandle&) = delete;
    FftwPlanHandle& operator=(const FftwPlanHandle&) = delete;
};

struct FftwPair {
    FftwApi& api;
    std::size_t n;
    FftwBuffer data;
    FftwPlanHandle forward;
    FftwPlanHandle inverse;
    FftwPair(FftwApi& owner, std::size_t size, unsigned flags)
        : api(owner), n(size), data(owner, 2 * size),
          forward(owner, owner.plan_dft(static_cast<int>(size), data.data, data.data, fftw_forward, flags)),
          inverse(owner, owner.plan_dft(static_cast<int>(size), data.data, data.data, fftw_backward, flags)) {}
    void initialize(const fftlab::Vector& input) {
        for (std::size_t i = 0; i < n; ++i) {
            data.data[2 * i] = input[i].real();
            data.data[2 * i + 1] = input[i].imag();
        }
    }
    void pair() {
        api.execute_dft(forward.plan, data.data, data.data);
        api.execute_dft(inverse.plan, data.data, data.data);
        const double scale = 1.0 / static_cast<double>(n);
        for (std::size_t i = 0; i < 2 * n; ++i) data.data[i] *= scale;
    }
};

volatile double sink = 0.0;

std::size_t calibrate(auto&& operation, double target_ms) {
    std::size_t iterations = 1;
    while (iterations < (1U << 20)) {
        const auto start = Clock::now();
        for (std::size_t i = 0; i < iterations; ++i) operation();
        const double elapsed = std::chrono::duration<double, std::milli>(Clock::now() - start).count();
        if (elapsed >= target_ms) break;
        iterations = std::min<std::size_t>(
            1U << 20,
            std::max(iterations + 1,
                     static_cast<std::size_t>(static_cast<double>(iterations) *
                                              std::clamp(target_ms / std::max(elapsed, 1e-9), 2.0, 16.0))));
    }
    return iterations;
}

struct Options {
    std::size_t n = 509;
    std::size_t samples = 31;
    std::size_t setup_samples = 1;
    std::size_t warmups = 5;
    double target_ms = 2.0;
    std::uint64_t seed = 20260812;
    bool info = false;
    bool self_test = false;
    bool raw_csv = false;
};

Options parse_options(int argc, char** argv) {
    Options options;
    auto value = [&](int& index) {
        if (++index >= argc) throw std::invalid_argument("missing option value");
        return std::string_view(argv[index]);
    };
    for (int i = 1; i < argc; ++i) {
        const std::string_view argument = argv[i];
        if (argument == "--info") options.info = true;
        else if (argument == "--self-test") options.self_test = true;
        else if (argument == "--raw-csv") options.raw_csv = true;
        else if (argument == "--size") options.n = std::stoull(std::string(value(i)));
        else if (argument == "--samples") options.samples = std::stoull(std::string(value(i)));
        else if (argument == "--setup-samples") options.setup_samples = std::stoull(std::string(value(i)));
        else if (argument == "--warmups") options.warmups = std::stoull(std::string(value(i)));
        else if (argument == "--target-ms") options.target_ms = std::stod(std::string(value(i)));
        else if (argument == "--seed") options.seed = std::stoull(std::string(value(i)));
        else throw std::invalid_argument("unknown option");
    }
    return options;
}

void emit(std::string_view phase, std::string_view backend, std::string_view algorithm,
          std::size_t n, std::size_t sample, std::size_t mode_order,
          std::size_t iterations, double ns_per_transform, std::size_t convolution_size,
          bool direct_cyclic) {
    std::cout << phase << ',' << backend << ',' << algorithm << ',' << n << ',' << sample << ','
              << mode_order << ',' << iterations << ',' << std::setprecision(12) << ns_per_transform
              << ',' << convolution_size << ',' << (direct_cyclic ? 1 : 0) << '\n';
}

void fftw_cross_check(FftwApi& api) {
    std::size_t checks = 0;
    std::mt19937_64 rng(0xF17A8B1ULL);
    std::uniform_real_distribution<double> dist(-1.0, 1.0);
    for (std::size_t n : {17U, 127U, 509U}) {
        fftlab::Vector input(n);
        for (auto& z : input) z = {dist(rng), dist(rng)};
        fftlab::RaderPlan reference_plan(n);
        fftlab::Vector reference(n), scratch(reference_plan.scratch_size());
        reference_plan.forward(input, reference, scratch);
        FftwPair fftw(api, n, fftw_estimate);
        fftw.initialize(input);
        api.execute_dft(fftw.forward.plan, fftw.data.data, fftw.data.data);
        for (std::size_t k = 0; k < n; ++k) {
            const fftlab::Complex got{fftw.data.data[2 * k], fftw.data.data[2 * k + 1]};
            ++checks;
            if (std::abs(got - reference[k]) > 5e-10 * (1.0 + std::abs(reference[k])))
                throw std::runtime_error("FFTW arbitrary-length cross-check mismatch");
        }
    }
    std::cout << "PASS: " << checks << " FFTW arbitrary-length cross-checks\n";
}

void benchmark(FftwApi& api, const Options& options) {
    if (options.n < 3 || !fftlab::is_prime(options.n))
        throw std::invalid_argument("arbitrary benchmark requires prime N >= 3");
    if (options.samples < 5 || options.setup_samples < 1 || options.target_ms <= 0)
        throw std::invalid_argument("bad benchmark parameters");

    std::mt19937_64 rng(options.seed ^ options.n);
    std::uniform_real_distribution<double> dist(-0.5, 0.5);
    fftlab::Vector input(options.n);
    for (auto& z : input) z = {dist(rng), dist(rng)};

    fftlab::BluesteinPlan blue_plan(options.n);
    fftlab::RaderPlan rader_plan(options.n);
    fftlab::Vector blue_a = input, blue_b(options.n), blue_scratch(blue_plan.scratch_size());
    fftlab::Vector rader_a = input, rader_b(options.n), rader_scratch(rader_plan.scratch_size());
    fftlab::Vector legacy_blue = input, legacy_rader = input;

    if (api.forget_wisdom) api.forget_wisdom();
    FftwPair estimate(api, options.n, fftw_estimate);
    if (api.forget_wisdom) api.forget_wisdom();
    FftwPair measure(api, options.n, fftw_measure);
    estimate.initialize(input);
    measure.initialize(input);

    auto pair_blue = [&] {
        blue_plan.forward(blue_a, blue_b, blue_scratch);
        blue_plan.inverse(blue_b, blue_a, blue_scratch);
    };
    auto pair_rader = [&] {
        rader_plan.forward(rader_a, rader_b, rader_scratch);
        rader_plan.inverse(rader_b, rader_a, rader_scratch);
    };
    auto pair_legacy_blue = [&] {
        legacy_blue = fftlab::bluestein(legacy_blue);
        legacy_blue = fftlab::bluestein(legacy_blue, true);
    };
    auto pair_legacy_rader = [&] {
        legacy_rader = fftlab::rader(legacy_rader);
        legacy_rader = fftlab::rader(legacy_rader, true);
    };

    for (std::size_t w = 0; w < options.warmups; ++w) {
        pair_blue(); pair_rader(); pair_legacy_blue(); pair_legacy_rader(); estimate.pair(); measure.pair();
    }
    const auto iterations = calibrate([&] { measure.pair(); }, options.target_ms);
    std::array<int, 6> modes{0, 1, 2, 3, 4, 5};

    if (options.raw_csv)
        std::cout << "phase,backend,algorithm,N,sample,mode_order,iterations,ns_per_transform,convolution_size,direct_cyclic\n";

    for (std::size_t sample = 0; sample < options.samples; ++sample) {
        std::shuffle(modes.begin(), modes.end(), rng);
        for (std::size_t order = 0; order < modes.size(); ++order) {
            const int mode = modes[order];
            const auto start = Clock::now();
            if (mode == 0) for (std::size_t i = 0; i < iterations; ++i) pair_legacy_blue();
            else if (mode == 1) for (std::size_t i = 0; i < iterations; ++i) pair_blue();
            else if (mode == 2) for (std::size_t i = 0; i < iterations; ++i) pair_legacy_rader();
            else if (mode == 3) for (std::size_t i = 0; i < iterations; ++i) pair_rader();
            else if (mode == 4) for (std::size_t i = 0; i < iterations; ++i) estimate.pair();
            else for (std::size_t i = 0; i < iterations; ++i) measure.pair();
            const double ns = std::chrono::duration<double, std::nano>(Clock::now() - start).count() /
                              static_cast<double>(2 * iterations);
            if (mode == 0) { sink = legacy_blue[0].real(); emit("execution", "fftlab", "legacy-bluestein", options.n, sample, order, iterations, ns, blue_plan.convolution_size(), false); }
            else if (mode == 1) { sink = blue_a[0].real(); emit("execution", "fftlab", "planned-bluestein", options.n, sample, order, iterations, ns, blue_plan.convolution_size(), false); }
            else if (mode == 2) { sink = legacy_rader[0].real(); emit("execution", "fftlab", "legacy-rader", options.n, sample, order, iterations, ns, rader_plan.convolution_size(), rader_plan.direct_cyclic_fft()); }
            else if (mode == 3) { sink = rader_a[0].real(); emit("execution", "fftlab", "planned-rader", options.n, sample, order, iterations, ns, rader_plan.convolution_size(), rader_plan.direct_cyclic_fft()); }
            else if (mode == 4) { sink = estimate.data.data[0]; emit("execution", "fftw", "estimate", options.n, sample, order, iterations, ns, 0, false); }
            else { sink = measure.data.data[0]; emit("execution", "fftw", "measure", options.n, sample, order, iterations, ns, 0, false); }
        }
    }

    FftwBuffer setup_buffer(api, 2 * options.n);
    std::array<int, 4> setup_modes{0, 1, 2, 3};
    for (std::size_t sample = 0; sample < options.setup_samples; ++sample) {
        std::shuffle(setup_modes.begin(), setup_modes.end(), rng);
        for (std::size_t order = 0; order < setup_modes.size(); ++order) {
            const int mode = setup_modes[order];
            if (mode >= 2 && api.forget_wisdom) api.forget_wisdom();
            const auto start = Clock::now();
            if (mode == 0) {
                fftlab::BluesteinPlan plan(options.n); sink = static_cast<double>(plan.convolution_size());
            } else if (mode == 1) {
                fftlab::RaderPlan plan(options.n); sink = static_cast<double>(plan.convolution_size());
            } else {
                const unsigned flags = mode == 2 ? fftw_estimate : fftw_measure;
                FftwPlanHandle fwd(api, api.plan_dft(static_cast<int>(options.n), setup_buffer.data, setup_buffer.data, fftw_forward, flags));
                FftwPlanHandle inv(api, api.plan_dft(static_cast<int>(options.n), setup_buffer.data, setup_buffer.data, fftw_backward, flags));
                sink = 1.0;
            }
            const double ns = std::chrono::duration<double, std::nano>(Clock::now() - start).count();
            if (mode == 0) emit("setup", "fftlab", "planned-bluestein", options.n, sample, order, 1, ns, blue_plan.convolution_size(), false);
            else if (mode == 1) emit("setup", "fftlab", "planned-rader", options.n, sample, order, 1, ns, rader_plan.convolution_size(), rader_plan.direct_cyclic_fft());
            else if (mode == 2) emit("setup", "fftw", "estimate", options.n, sample, order, 1, ns, 0, false);
            else emit("setup", "fftw", "measure", options.n, sample, order, 1, ns, 0, false);
        }
    }
}

} // namespace

int main(int argc, char** argv) {
    try {
        const auto options = parse_options(argc, argv);
        FftwApi api;
        if (options.info) {
            std::cout << "fftw_available=" << (api ? "yes" : "no") << '\n';
            std::cout << "fftw_library=" << api.library << '\n';
            std::cout << "fftw_version=" << api.version << '\n';
            return api ? 0 : 77;
        }
        if (options.self_test) {
            fftlab::arbitrary_plan_tests();
            if (api) fftw_cross_check(api);
            else std::cout << "SKIP: FFTW runtime unavailable\n";
            return 0;
        }
        if (!api) { std::cerr << "FFTW runtime not found\n"; return 77; }
        if (options.raw_csv) { benchmark(api, options); return 0; }
        std::cerr << "use --info, --self-test, or --raw-csv\n";
        return 2;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return 2;
    }
}
