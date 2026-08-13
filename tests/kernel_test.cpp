#include "fftlab/kernel.hpp"
#include "fftlab/plan.hpp"
#include <iostream>
#include <random>
using namespace fftlab;
int main(){
  auto caps=kernel_capabilities(); std::size_t checks=0;
  std::mt19937_64 rng(7); std::uniform_real_distribution<double>d(-1,1);
  for(std::size_t n:{1u,2u,4u,8u,16u,64u,256u,1024u}){
    Vector64 x(n); for(auto&z:x)z={d(rng),d(rng)};
    Radix2Plan<double> refp(n); auto ref=x; refp.forward_inplace(ref);
    std::vector<KernelIsa> modes{KernelIsa::Scalar}; if(caps.avx2)modes.push_back(KernelIsa::Avx2); if(caps.avx512)modes.push_back(KernelIsa::Avx512);
    for(auto isa:modes){KernelRadix2Plan p(n,isa);auto y=x;p.forward_inplace(y);for(std::size_t i=0;i<n;++i){++checks;if(std::abs(y[i]-ref[i])>3e-12*(1+std::abs(ref[i])))throw std::runtime_error("kernel mismatch");}p.inverse_inplace(y);for(std::size_t i=0;i<n;++i){++checks;if(std::abs(y[i]-x[i])>3e-12*(1+std::abs(x[i])))throw std::runtime_error("kernel rt");}}
  }
  std::cout<<"PASS: "<<checks<<" kernel checks; avx2="<<caps.avx2<<" avx512="<<caps.avx512<<"\n";
}
