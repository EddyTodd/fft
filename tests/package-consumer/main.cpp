#include <fftlab/fftlab.hpp>
#include <fftlab/version.hpp>

int main() {
    return fftlab::version_major == 1 ? 0 : 1;
}
