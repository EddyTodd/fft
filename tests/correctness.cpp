#include "fftlab/planner.hpp"
#include "fftlab/oracle.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <numbers>
#include <random>
#include <string_view>

using namespace fftlab;

static std::size_t checks = 0;
void require(bool ok, std::string_view msg) { ++checks; if (!ok) throw std::runtime_error(std::string(msg)); }

template<FftScalar T>
T tolerance() { return std::same_as<T,float> ? static_cast<T>(8e-4) : static_cast<T>(8e-10); }

template<FftScalar T>
bool approx(ComplexT<T> a, ComplexT<T> b, T mult=T{1}) {
    return std::abs(a-b) <= tolerance<T>()*mult*(T{1}+std::abs(b));
}
template<FftScalar T>
bool approx_real(T a, T b, T mult=T{1}) {
    return std::abs(a-b) <= tolerance<T>()*mult*(T{1}+std::abs(b));
}

template<FftScalar T>
void oracle_matrix() {
    std::mt19937_64 rng(0xF17F17);
    std::uniform_real_distribution<double> dist(-1,1);
    for (std::size_t n : {0u,1u,2u,3u,4u,5u,6u,7u,8u,9u,10u,12u,15u,16u,17u,18u,20u,21u,25u,27u,31u,35u,49u,60u,64u,75u,121u,127u}) {
        VectorT<T> x(n); for(auto& z:x) z={static_cast<T>(dist(rng)),static_cast<T>(dist(rng))};
        auto ref=oracle_dft<T>(x);
        auto auto_y=transform(x);
        for(std::size_t i=0;i<n;++i) require(std::abs(OracleComplex{auto_y[i].real(),auto_y[i].imag()}-ref[i]) <= static_cast<long double>(tolerance<T>()*T{4})*(1.0L+std::abs(ref[i])),"auto algorithm oracle");
        Plan<T> plan(n); VectorT<T> y(n), scratch(plan.scratch_size()); plan.forward(x,y,scratch);
        for(std::size_t i=0;i<n;++i) require(std::abs(OracleComplex{y[i].real(),y[i].imag()}-ref[i]) <= static_cast<long double>(tolerance<T>()*T{4})*(1.0L+std::abs(ref[i])),"structural plan oracle");
        VectorT<T> back(n); plan.inverse(y,back,scratch);
        for(std::size_t i=0;i<n;++i) require(approx(back[i],x[i],T{3}),"plan roundtrip");
    }
}

template<FftScalar T>
void catalog_oracle() {
    std::mt19937_64 rng(0xCA7A106);
    std::uniform_real_distribution<double> dist(-1,1);
    for (std::size_t n : {0u,1u,2u,3u,4u,5u,6u,7u,8u,10u,12u,15u,16u,17u,20u,21u,25u,31u,32u}) {
        VectorT<T> x(n);
        for (auto& z : x) z={static_cast<T>(dist(rng)),static_cast<T>(dist(rng))};
        const auto ref = oracle_dft<T>(x);
        for (const auto algo : all_algos) {
            if (algo == Algo::Auto || algo == Algo::Dft || !supports(algo,n)) continue;
            const auto y = transform(x, algo, Direction::Forward);
            require(y.size()==n,"catalog output size");
            for (std::size_t i=0;i<n;++i) {
                const OracleComplex got{static_cast<long double>(y[i].real()),static_cast<long double>(y[i].imag())};
                require(std::abs(got-ref[i]) <= static_cast<long double>(tolerance<T>()*T{10})*(1.0L+std::abs(ref[i])),"catalog forward oracle");
            }
            const auto back = transform(y, algo, Direction::Inverse);
            for (std::size_t i=0;i<n;++i) require(approx(back[i],x[i],T{8}),"catalog inverse roundtrip");
        }
    }
}

template<FftScalar T>
void codelet_tests() {
    for (const std::size_t radix : {2u,3u,4u,5u,7u}) {
        SmallDftCodelet<T> codelet(radix);
        require(codelet.radix()==radix,"codelet radix");
        if(radix==2 || radix==4) require(codelet.stored_roots()==0,"specialized codelet has no root matrix");
        else require(codelet.stored_roots()==radix*radix,"planned codelet root matrix");
        VectorT<T> x(radix);
        for(std::size_t i=0;i<radix;++i) x[i]={static_cast<T>(i+1)/static_cast<T>(radix),static_cast<T>(2*i+1)/static_cast<T>(radix+1)};
        const auto ref=oracle_dft<T>(x);
        auto y=x; codelet.execute(y,Direction::Forward);
        for(std::size_t i=0;i<radix;++i) {
            const OracleComplex got{static_cast<long double>(y[i].real()),static_cast<long double>(y[i].imag())};
            require(std::abs(got-ref[i]) <= static_cast<long double>(tolerance<T>()*T{4})*(1.0L+std::abs(ref[i])),"codelet oracle");
        }
        codelet.execute(y,Direction::Inverse);
        for(std::size_t i=0;i<radix;++i) require(approx(y[i],x[i],T{4}),"codelet roundtrip");
    }
}

template<FftScalar T>
void identities() {
    for(std::size_t n:{3u,8u,12u,15u,17u,25u,31u,64u}) {
        VectorT<T> impulse(n); impulse[0]={1,0}; auto f=transform(impulse);
        for(auto z:f) require(approx(z,{1,0},T{2}),"impulse identity");
        VectorT<T> constant(n,{T{1},T{0}}); f=transform(constant);
        require(approx(f[0],{static_cast<T>(n),0},T{2}),"constant DC");
        for(std::size_t k=1;k<n;++k) require(std::abs(f[k]) <= tolerance<T>()*T{8}*static_cast<T>(n),"constant non-DC");
        const std::size_t tone= n>1 ? std::min<std::size_t>(2,n-1) : 0;
        VectorT<T> single(n);
        for(std::size_t j=0;j<n;++j) {
            const T angle=T{2}*std::numbers::pi_v<T>*static_cast<T>(tone*j)/static_cast<T>(n);
            single[j]=root<T>(angle);
        }
        f=transform(single);
        for(std::size_t k=0;k<n;++k) {
            if(k==tone) require(approx(f[k],{static_cast<T>(n),0},T{3}),"single tone bin");
            else require(std::abs(f[k]) <= tolerance<T>()*T{12}*static_cast<T>(n),"single tone leakage");
        }
    }
}

template<FftScalar T>
void real_tests() {
    std::mt19937_64 rng(2026); std::uniform_real_distribution<double> dist(-1,1);
    for(std::size_t n:{1u,2u,4u,8u,16u,64u}) {
        RealRadix2Plan<T> plan(n);
        std::vector<T> x(n); for(auto& v:x) v=static_cast<T>(dist(rng));
        std::vector<ComplexT<T>> spectrum(plan.spectrum_size()), scratch(plan.scratch_size());
        plan.forward(x,spectrum,scratch);
        VectorT<T> complex_x(n); for(std::size_t i=0;i<n;++i) complex_x[i]={x[i],0};
        auto full=dft(complex_x);
        for(std::size_t k=0;k<spectrum.size();++k) require(approx(spectrum[k],full[k],T{3}),"real half spectrum");
        if(n>1) {
            for(std::size_t k=1;k<n/2;++k) require(approx(full[n-k],std::conj(full[k]),T{3}),"Hermitian symmetry");
        }
        std::vector<T> back(n); plan.inverse(spectrum,back,scratch);
        for(std::size_t i=0;i<n;++i) require(approx_real(back[i],x[i],T{3}),"real roundtrip");
    }
}

template<FftScalar T>
void explicit_mechanisms() {
    std::mt19937_64 rng(99); std::uniform_real_distribution<double> dist(-1,1);
    for(std::size_t n:{12u,18u,20u,25u,27u,45u,49u,60u,75u,121u}) {
        VectorT<T> x(n); for(auto& z:x)z={static_cast<T>(dist(rng)),static_cast<T>(dist(rng))}; auto ref=oracle_dft<T>(x);
        MixedRadixPlan<T> p(n); VectorT<T> y=x,s(p.scratch_size()); p.forward_inplace(y,s);
        for(std::size_t i=0;i<n;++i) require(std::abs(OracleComplex{y[i].real(),y[i].imag()}-ref[i]) <= static_cast<long double>(tolerance<T>()*T{6})*(1.0L+std::abs(ref[i])),"planned mixed radix");
    }
    for(std::size_t n:{6u,10u,12u,15u,18u,20u,21u,35u,45u,60u,75u,143u}) {
        auto split=coprime_factor_split(n); if(split.first==0) continue;
        VectorT<T> x(n); for(auto& z:x)z={static_cast<T>(dist(rng)),static_cast<T>(dist(rng))}; auto ref=oracle_dft<T>(x);
        GoodThomasPlan<T> p(n); require(p.twiddle_count()==0,"PFA top-level is twiddle free");
        VectorT<T> y=x,s(p.scratch_size()); p.forward_inplace(y,s);
        for(std::size_t i=0;i<n;++i) require(std::abs(OracleComplex{y[i].real(),y[i].imag()}-ref[i]) <= static_cast<long double>(tolerance<T>()*T{8})*(1.0L+std::abs(ref[i])),"Good Thomas plan");
    }
    for(std::size_t n:{8u,16u,32u,64u,128u}) {
        VectorT<T> x(n); for(auto& z:x)z={static_cast<T>(dist(rng)),static_cast<T>(dist(rng))};
        auto ref=oracle_dft<T>(x); auto y=modified_split_radix(x);
        for(std::size_t i=0;i<n;++i) require(std::abs(OracleComplex{y[i].real(),y[i].imag()}-ref[i]) <= static_cast<long double>(tolerance<T>()*T{8})*(1.0L+std::abs(ref[i])),"modified split radix");
    }
}

void planner_policy() {
    require(Plan<double>(0).algorithm()==PlanAlgorithm::Identity,"N0 identity");
    require(Plan<double>(1).algorithm()==PlanAlgorithm::Identity,"N1 identity");
    require(Plan<double>(64).algorithm()==PlanAlgorithm::Radix2,"power2 radix2");
    require(Plan<double>(25).algorithm()==PlanAlgorithm::MixedRadix,"prime power mixed");
    require(Plan<double>(15).algorithm()==PlanAlgorithm::GoodThomas,"coprime PFA");
    require(Plan<double>(17).algorithm()==PlanAlgorithm::Rader,"short Rader");
    require(Plan<double>(31).algorithm()==PlanAlgorithm::Bluestein,"equal-convolution prime Bluestein");
    require(Plan<double>(121).algorithm()==PlanAlgorithm::Bluestein,"rough composite Bluestein fallback");
    bool threw=false; try{Plan<double> p(12,{PlanPreference::Rader,true});}catch(const std::invalid_argument&){threw=true;} require(threw,"forced invalid policy rejects");
}

int main(){
    oracle_matrix<float>(); oracle_matrix<double>();
    catalog_oracle<float>(); catalog_oracle<double>();
    codelet_tests<float>(); codelet_tests<double>();
    identities<float>(); identities<double>();
    real_tests<float>(); real_tests<double>();
    explicit_mechanisms<float>(); explicit_mechanisms<double>();
    planner_policy();
    std::cout<<"PASS: "<<checks<<" checks\n";
}
