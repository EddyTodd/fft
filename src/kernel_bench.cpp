#include "fftlab/kernel.hpp"
#include "fftlab/plan.hpp"

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
void* load_symbol(LibraryHandle handle, const char* name) {
    return reinterpret_cast<void*>(GetProcAddress(handle, name));
}
void close_library(LibraryHandle handle) {
    if (handle) FreeLibrary(handle);
}
#else
using LibraryHandle = void*;
LibraryHandle open_library(const char* name) { return dlopen(name, RTLD_NOW | RTLD_LOCAL); }
void* load_symbol(LibraryHandle handle, const char* name) { return dlsym(handle, name); }
void close_library(LibraryHandle handle) {
    if (handle) dlclose(handle);
}
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
            if (handle) {
                library = candidate;
                break;
            }
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

    ~FftwApi() {
        if (cleanup) cleanup();
        close_library(handle);
    }
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

struct FftwComplexPair {
    FftwApi& api;
    std::size_t n;
    FftwBuffer data;
    FftwPlanHandle forward;
    FftwPlanHandle inverse;

    FftwComplexPair(FftwApi& owner, std::size_t size, unsigned flags)
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
    std::size_t n = 1024;
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

void emit(std::string_view phase, std::string_view backend, std::string_view policy,
          std::size_t n, std::size_t sample, std::size_t mode_order,
          std::size_t iterations, double ns_per_transform) {
    std::cout << phase << ',' << backend << ',' << policy << ',' << n << ',' << sample << ','
              << mode_order << ',' << iterations << ',' << std::setprecision(12)
              << ns_per_transform << '\n';
}

void fftw_cross_check(FftwApi& api) {
    std::size_t checks = 0;
    for (std::size_t n : {8U, 64U, 256U}) {
        fftlab::Vector input(n);
        for (std::size_t i = 0; i < n; ++i) {
            input[i] = {static_cast<double>(i % 17) / 17.0 - 0.5,
                        static_cast<double>(i % 19) / 19.0 - 0.5};
        }

        fftlab::KernelRadix2Plan reference_plan(n, fftlab::KernelIsa::Scalar);
        auto reference = input;
        reference_plan.forward_inplace(reference);

        FftwComplexPair fftw(api, n, fftw_estimate);
        fftw.initialize(input);
        api.execute_dft(fftw.forward.plan, fftw.data.data, fftw.data.data);
        for (std::size_t i = 0; i < n; ++i) {
            ++checks;
            const fftlab::Complex actual{fftw.data.data[2 * i], fftw.data.data[2 * i + 1]};
            if (std::abs(actual - reference[i]) > 2e-10 * (1.0 + std::abs(reference[i]))) {
                throw std::runtime_error("FFTW/kernel frequency-domain mismatch");
            }
        }
    }
    std::cout << "PASS: " << checks << " FFTW kernel cross-checks\n";
}

void benchmark(FftwApi& api, const Options& options) {
    if (options.n < 2 || !fftlab::pow2(options.n))
        throw std::invalid_argument("kernel benchmark requires power-of-two N >= 2");
    if (options.samples < 5 || options.setup_samples < 1 || options.target_ms <= 0.0)
        throw std::invalid_argument("bad benchmark parameters");

    const auto capabilities = fftlab::kernel_capabilities();
    if (!capabilities.avx2 || !capabilities.avx512)
        throw std::runtime_error("formal kernel matrix requires AVX2/FMA and AVX-512/FMA");

    std::mt19937_64 rng(options.seed ^ options.n);
    std::uniform_real_distribution<double> distribution(-0.5, 0.5);
    fftlab::Vector input(options.n);
    for (auto& value : input) value = {distribution(rng), distribution(rng)};

    fftlab::Radix2Plan baseline(options.n);
    fftlab::KernelRadix2Plan scalar(options.n, fftlab::KernelIsa::Scalar);
    fftlab::KernelRadix2Plan avx2(options.n, fftlab::KernelIsa::Avx2);
    fftlab::KernelRadix2Plan avx512(options.n, fftlab::KernelIsa::Avx512);
    fftlab::KernelRadix2Plan automatic(options.n, fftlab::KernelIsa::Auto);

    if (api.forget_wisdom) api.forget_wisdom();
    FftwComplexPair estimate(api, options.n, fftw_estimate);
    if (api.forget_wisdom) api.forget_wisdom();
    FftwComplexPair measure(api, options.n, fftw_measure);

    auto baseline_buffer = input;
    auto scalar_buffer = input;
    auto avx2_buffer = input;
    auto avx512_buffer = input;
    auto automatic_buffer = input;
    estimate.initialize(input);
    measure.initialize(input);

    auto pair = [](auto& plan, auto& buffer) {
        plan.forward_inplace(buffer);
        plan.inverse_inplace(buffer);
    };

    for (std::size_t warmup = 0; warmup < options.warmups; ++warmup) {
        pair(baseline, baseline_buffer);
        pair(scalar, scalar_buffer);
        pair(avx2, avx2_buffer);
        pair(avx512, avx512_buffer);
        pair(automatic, automatic_buffer);
        estimate.pair();
        measure.pair();
    }

    const auto iterations = calibrate([&] { pair(avx2, avx2_buffer); }, options.target_ms);
    std::array<int, 7> modes{0, 1, 2, 3, 4, 5, 6};
    std::cout << "phase,backend,policy,N,sample,mode_order,iterations,ns_per_transform\n";

    for (std::size_t sample = 0; sample < options.samples; ++sample) {
        std::shuffle(modes.begin(), modes.end(), rng);
        for (std::size_t order = 0; order < modes.size(); ++order) {
            const int mode = modes[order];
            const auto start = Clock::now();
            if (mode == 0) {
                for (std::size_t i = 0; i < iterations; ++i) pair(baseline, baseline_buffer);
            } else if (mode == 1) {
                for (std::size_t i = 0; i < iterations; ++i) pair(scalar, scalar_buffer);
            } else if (mode == 2) {
                for (std::size_t i = 0; i < iterations; ++i) pair(avx2, avx2_buffer);
            } else if (mode == 3) {
                for (std::size_t i = 0; i < iterations; ++i) pair(avx512, avx512_buffer);
            } else if (mode == 4) {
                for (std::size_t i = 0; i < iterations; ++i) pair(automatic, automatic_buffer);
            } else if (mode == 5) {
                for (std::size_t i = 0; i < iterations; ++i) estimate.pair();
            } else {
                for (std::size_t i = 0; i < iterations; ++i) measure.pair();
            }
            const double ns = std::chrono::duration<double, std::nano>(Clock::now() - start).count() /
                              static_cast<double>(2 * iterations);

            if (mode == 0) {
                sink = baseline_buffer[0].real();
                emit("execution", "fftlab-plan", "legacy", options.n, sample, order, iterations, ns);
            } else if (mode == 1) {
                sink = scalar_buffer[0].real();
                emit("execution", "kernel", "scalar", options.n, sample, order, iterations, ns);
            } else if (mode == 2) {
                sink = avx2_buffer[0].real();
                emit("execution", "kernel", "avx2", options.n, sample, order, iterations, ns);
            } else if (mode == 3) {
                sink = avx512_buffer[0].real();
                emit("execution", "kernel", "avx512", options.n, sample, order, iterations, ns);
            } else if (mode == 4) {
                sink = automatic_buffer[0].real();
                const std::string policy = "auto->" + std::string(fftlab::kernel_name(automatic.selected_isa()));
                emit("execution", "kernel", policy, options.n, sample, order, iterations, ns);
            } else if (mode == 5) {
                sink = estimate.data.data[0];
                emit("execution", "fftw", "estimate", options.n, sample, order, iterations, ns);
            } else {
                sink = measure.data.data[0];
                emit("execution", "fftw", "measure", options.n, sample, order, iterations, ns);
            }
        }
    }

    std::array<int, 7> setup_modes{0, 1, 2, 3, 4, 5, 6};
    FftwBuffer setup_buffer(api, 2 * options.n);
    for (std::size_t sample = 0; sample < options.setup_samples; ++sample) {
        std::shuffle(setup_modes.begin(), setup_modes.end(), rng);
        for (std::size_t order = 0; order < setup_modes.size(); ++order) {
            const int mode = setup_modes[order];
            if (mode >= 5 && api.forget_wisdom) api.forget_wisdom();
            const auto start = Clock::now();
            std::string policy;
            if (mode == 0) {
                fftlab::Radix2Plan plan(options.n);
                sink = static_cast<double>(plan.stored_twiddles());
                policy = "legacy";
            } else if (mode == 1) {
                fftlab::KernelRadix2Plan plan(options.n, fftlab::KernelIsa::Scalar);
                sink = static_cast<double>(plan.stored_twiddles());
                policy = "scalar";
            } else if (mode == 2) {
                fftlab::KernelRadix2Plan plan(options.n, fftlab::KernelIsa::Avx2);
                sink = static_cast<double>(plan.stored_twiddles());
                policy = "avx2";
            } else if (mode == 3) {
                fftlab::KernelRadix2Plan plan(options.n, fftlab::KernelIsa::Avx512);
                sink = static_cast<double>(plan.stored_twiddles());
                policy = "avx512";
            } else if (mode == 4) {
                fftlab::KernelRadix2Plan plan(options.n, fftlab::KernelIsa::Auto);
                sink = plan.tuning_ns();
                policy = "auto->" + std::string(fftlab::kernel_name(plan.selected_isa()));
            } else {
                const unsigned flags = mode == 5 ? fftw_estimate : fftw_measure;
                FftwPlanHandle forward(api, api.plan_dft(static_cast<int>(options.n), setup_buffer.data,
                                                         setup_buffer.data, fftw_forward, flags));
                FftwPlanHandle inverse(api, api.plan_dft(static_cast<int>(options.n), setup_buffer.data,
                                                         setup_buffer.data, fftw_backward, flags));
                sink = 1.0;
                policy = mode == 5 ? "estimate" : "measure";
            }
            const double ns = std::chrono::duration<double, std::nano>(Clock::now() - start).count();
            const std::string_view backend = mode == 0 ? "fftlab-plan" : (mode < 5 ? "kernel" : "fftw");
            emit("setup", backend, policy, options.n, sample, order, 1, ns);
        }
    }
}

} // namespace

int main(int argc, char** argv) {
    try {
        const auto options = parse_options(argc, argv);
        FftwApi fftw;
        if (options.info) {
            const auto capabilities = fftlab::kernel_capabilities();
            std::cout << "fftw_available=" << (fftw ? "yes" : "no") << '\n'
                      << "fftw_library=" << fftw.library << '\n'
                      << "fftw_version=" << fftw.version << '\n'
                      << "avx2_fma=" << (capabilities.avx2 ? "yes" : "no") << '\n'
                      << "avx512_fma=" << (capabilities.avx512 ? "yes" : "no") << '\n';
            return 0;
        }
        if (options.self_test) {
            fftlab::kernel_tests();
            if (fftw) fftw_cross_check(fftw);
            else std::cout << "SKIP: FFTW runtime not found; kernel correctness still passed\n";
            return 0;
        }
        if (options.raw_csv) {
            if (!fftw) {
                std::cerr << "FFTW runtime not found\n";
                return 77;
            }
            benchmark(fftw, options);
            return 0;
        }
        std::cerr << "use --info, --self-test, or --raw-csv\n";
        return 2;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return 2;
    }
}
