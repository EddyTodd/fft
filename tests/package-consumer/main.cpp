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

int main() {
    return fftlab::version_major == 1 ? 0 : 1;
}
