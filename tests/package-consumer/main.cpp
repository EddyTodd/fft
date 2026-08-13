#include <fftlab/arbitrary_algorithms.hpp>
#include <fftlab/arbitrary_plan.hpp>
#include <fftlab/codelet.hpp>
#include <fftlab/convolution_plan.hpp>
#include <fftlab/fft.hpp>
#include <fftlab/fftlab.hpp>
#include <fftlab/good_thomas_plan.hpp>
#include <fftlab/kernel.hpp>
#include <fftlab/mixed_plan.hpp>
#include <fftlab/oracle.hpp>
#include <fftlab/plan.hpp>
#include <fftlab/planner.hpp>
#include <fftlab/power2_algorithms.hpp>
#include <fftlab/types.hpp>
#include <fftlab/version.hpp>

#include <array>
#include <cmath>
#include <span>

int main() {
    if (fftlab::version_major != 1) return 1;

    std::array<fftlab::Complex64, 4> data{
        fftlab::Complex64{1.0, 0.0},
        fftlab::Complex64{0.0, 0.0},
        fftlab::Complex64{0.0, 0.0},
        fftlab::Complex64{0.0, 0.0},
    };
    fftlab::KernelRadix2Plan plan(4, fftlab::KernelIsa::Scalar);
    plan.forward_inplace(std::span<fftlab::Complex64>{data});
    for (const auto& value : data) {
        if (std::abs(value.real() - 1.0) > 1e-12 || std::abs(value.imag()) > 1e-12) return 2;
    }
    plan.inverse_inplace(std::span<fftlab::Complex64>{data});
    if (std::abs(data[0].real() - 1.0) > 1e-12 || std::abs(data[0].imag()) > 1e-12) return 3;
    for (std::size_t i = 1; i < data.size(); ++i) {
        if (std::abs(data[i]) > 1e-12) return 4;
    }
    return 0;
}
