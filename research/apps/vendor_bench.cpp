#include "fftlab/plan.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#ifdef _WIN32
#include <windows.h>
#else
#include <dlfcn.h>
#endif
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <vector>

namespace {
using Clock = std::chrono::steady_clock;
using FPlan = void*;
constexpr int FFTW_FORWARD = -1;
constexpr int FFTW_BACKWARD = 1;
constexpr unsigned FFTW_MEASURE = 0U;
constexpr unsigned FFTW_ESTIMATE = 1U << 6;

#ifdef _WIN32
using LibraryHandle = HMODULE;
LibraryHandle open_library(const char* name) { return LoadLibraryA(name); }
void* load_symbol(LibraryHandle h, const char* name) { return reinterpret_cast<void*>(GetProcAddress(h, name)); }
void close_library(LibraryHandle h) { if (h) FreeLibrary(h); }
#else
using LibraryHandle = void*;
LibraryHandle open_library(const char* name) { return dlopen(name, RTLD_NOW | RTLD_LOCAL); }
void* load_symbol(LibraryHandle h, const char* name) { return dlsym(h, name); }
void close_library(LibraryHandle h) { if (h) dlclose(h); }
#endif

struct FftwApi {
    using PlanDft = FPlan(*)(int, double*, double*, int, unsigned);
    using PlanR2C = FPlan(*)(int, double*, double*, unsigned);
    using PlanC2R = FPlan(*)(int, double*, double*, unsigned);
    using ExecuteDft = void(*)(FPlan, double*, double*);
    using ExecuteR2C = void(*)(FPlan, double*, double*);
    using ExecuteC2R = void(*)(FPlan, double*, double*);
    using Destroy = void(*)(FPlan);
    using Malloc = void*(*)(std::size_t);
    using Free = void(*)(void*);
    using ForgetWisdom = void(*)();
    using Cleanup = void(*)();

    LibraryHandle handle{};
    std::string library, version;
    PlanDft plan_dft{}; PlanR2C plan_r2c{}; PlanC2R plan_c2r{};
    ExecuteDft execute_dft{}; ExecuteR2C execute_r2c{}; ExecuteC2R execute_c2r{};
    Destroy destroy_plan{}; Malloc malloc_fn{}; Free free_fn{}; ForgetWisdom forget_wisdom{}; Cleanup cleanup{};

    FftwApi() {
        for (const char* candidate : {"libfftw3.so.3", "libfftw3.so", "libfftw3.3.dylib", "libfftw3.dylib", "libfftw3-3.dll", "fftw3.dll"}) {
            handle = open_library(candidate);
            if (handle) { library = candidate; break; }
        }
        if (!handle) return;
        auto sym = [&](auto& fn, const char* name) { fn = reinterpret_cast<std::remove_reference_t<decltype(fn)>>(load_symbol(handle, name)); if (!fn) throw std::runtime_error(std::string("missing FFTW symbol: ") + name); };
        sym(plan_dft, "fftw_plan_dft_1d"); sym(plan_r2c, "fftw_plan_dft_r2c_1d"); sym(plan_c2r, "fftw_plan_dft_c2r_1d");
        sym(execute_dft, "fftw_execute_dft"); sym(execute_r2c, "fftw_execute_dft_r2c"); sym(execute_c2r, "fftw_execute_dft_c2r");
        sym(destroy_plan, "fftw_destroy_plan"); sym(malloc_fn, "fftw_malloc"); sym(free_fn, "fftw_free");
        forget_wisdom = reinterpret_cast<ForgetWisdom>(load_symbol(handle, "fftw_forget_wisdom"));
        cleanup = reinterpret_cast<Cleanup>(load_symbol(handle, "fftw_cleanup"));
        if (const auto* v = static_cast<const char*>(load_symbol(handle, "fftw_version"))) version = v;
    }
    ~FftwApi() { if (cleanup) cleanup(); close_library(handle); }
    FftwApi(const FftwApi&) = delete; FftwApi& operator=(const FftwApi&) = delete;
    explicit operator bool() const noexcept { return handle != nullptr; }
};

struct FftwBuffer {
    FftwApi* api{}; double* p{};
    FftwBuffer() = default;
    FftwBuffer(FftwApi& a, std::size_t doubles) : api(&a), p(static_cast<double*>(a.malloc_fn(sizeof(double) * doubles))) { if (!p) throw std::bad_alloc(); }
    ~FftwBuffer() { if (p) api->free_fn(p); }
    FftwBuffer(const FftwBuffer&) = delete; FftwBuffer& operator=(const FftwBuffer&) = delete;
};
struct FftwPlanHandle {
    FftwApi* api{}; FPlan p{};
    FftwPlanHandle() = default; FftwPlanHandle(FftwApi& a, FPlan q) : api(&a), p(q) { if (!p) throw std::runtime_error("FFTW plan creation failed"); }
    ~FftwPlanHandle() { if (p) api->destroy_plan(p); }
    FftwPlanHandle(const FftwPlanHandle&) = delete; FftwPlanHandle& operator=(const FftwPlanHandle&) = delete;
};

struct FftwComplexPair {
    FftwApi& api; std::size_t n; FftwBuffer data; FftwPlanHandle fwd, inv;
    FftwComplexPair(FftwApi& a, std::size_t size, unsigned flags) : api(a), n(size), data(a, 2 * size),
        fwd(a, a.plan_dft(int(size), data.p, data.p, FFTW_FORWARD, flags)), inv(a, a.plan_dft(int(size), data.p, data.p, FFTW_BACKWARD, flags)) {}
    void init(const fftlab::Vector& x) { for (std::size_t i=0;i<n;++i){data.p[2*i]=x[i].real();data.p[2*i+1]=x[i].imag();} }
    void pair() { api.execute_dft(fwd.p, data.p, data.p); api.execute_dft(inv.p, data.p, data.p); const double s=1.0/double(n); for(std::size_t i=0;i<2*n;++i)data.p[i]*=s; }
};
struct FftwRealPair {
    FftwApi& api; std::size_t n; FftwBuffer real, spectrum; FftwPlanHandle fwd, inv;
    FftwRealPair(FftwApi& a, std::size_t size, unsigned flags) : api(a), n(size), real(a,size), spectrum(a,2*(size/2+1)),
        fwd(a,a.plan_r2c(int(size),real.p,spectrum.p,flags)), inv(a,a.plan_c2r(int(size),spectrum.p,real.p,flags)) {}
    void init(const fftlab::RealVector& x) { std::copy(x.begin(),x.end(),real.p); }
    void pair() { api.execute_r2c(fwd.p,real.p,spectrum.p); api.execute_c2r(inv.p,spectrum.p,real.p); const double s=1.0/double(n); for(std::size_t i=0;i<n;++i)real.p[i]*=s; }
};

volatile double sink = 0;
std::size_t calibrate(auto&& fn, double target_ms) {
    std::size_t it=1;
    while(it<(1u<<20)) { auto t=Clock::now(); for(std::size_t i=0;i<it;++i) fn(); double ms=std::chrono::duration<double,std::milli>(Clock::now()-t).count(); if(ms>=target_ms) break; it=std::min<std::size_t>(1u<<20,std::max(it+1,std::size_t(double(it)*std::clamp(target_ms/std::max(ms,1e-9),2.0,16.0)))); }
    return it;
}

struct Opt { std::size_t n=1024,samples=31,setup_samples=3,warmups=5; double target_ms=2; std::uint64_t seed=20260812; bool info=false,self=false,raw=false; };
Opt parse(int ac,char**av){ Opt o; auto val=[&](int&i){if(++i>=ac)throw std::invalid_argument("missing value");return std::string_view(av[i]);}; for(int i=1;i<ac;++i){std::string_view a=av[i]; if(a=="--info")o.info=true; else if(a=="--self-test")o.self=true; else if(a=="--raw-csv")o.raw=true; else if(a=="--size")o.n=std::stoull(std::string(val(i))); else if(a=="--samples")o.samples=std::stoull(std::string(val(i))); else if(a=="--setup-samples")o.setup_samples=std::stoull(std::string(val(i))); else if(a=="--warmups")o.warmups=std::stoull(std::string(val(i))); else if(a=="--target-ms")o.target_ms=std::stod(std::string(val(i))); else if(a=="--seed")o.seed=std::stoull(std::string(val(i))); else throw std::invalid_argument("unknown option");} return o; }

void self_test(FftwApi& api) {
    std::mt19937_64 rng(0xFF7AULL); std::uniform_real_distribution<double>d(-1,1); std::size_t checks=0; auto req=[&](bool x,const char* m){++checks;if(!x)throw std::runtime_error(m);};
    for(std::size_t n:{8u,64u,256u}){
        fftlab::Vector x(n);for(auto&z:x)z={d(rng),d(rng)}; fftlab::Radix2Plan own(n); auto ref=x; own.forward_inplace(ref);
        FftwComplexPair fp(api,n,FFTW_ESTIMATE); fp.init(x); api.execute_dft(fp.fwd.p,fp.data.p,fp.data.p);
        for(std::size_t i=0;i<n;++i) req(std::abs(fftlab::Complex{fp.data.p[2*i],fp.data.p[2*i+1]}-ref[i])<2e-10*(1+std::abs(ref[i])),"FFTW complex mismatch");
        fftlab::RealVector r(n);for(auto&v:r)v=d(rng); fftlab::RealRadix2Plan rp(n); fftlab::HalfSpectrum hs(n/2+1); fftlab::Vector scratch(n/2); rp.forward(r,hs,scratch);
        FftwRealPair fr(api,n,FFTW_ESTIMATE); fr.init(r); api.execute_r2c(fr.fwd.p,fr.real.p,fr.spectrum.p);
        for(std::size_t k=0;k<=n/2;++k) req(std::abs(fftlab::Complex{fr.spectrum.p[2*k],fr.spectrum.p[2*k+1]}-hs[k])<2e-10*(1+std::abs(hs[k])),"FFTW real mismatch");
    }
    std::cout<<"PASS: "<<checks<<" FFTW cross-checks\n";
}

void emit(std::string_view phase,std::string_view backend,std::string_view planner,std::string_view kind,std::size_t n,std::size_t sample,std::size_t it,double ns){std::cout<<phase<<','<<backend<<','<<planner<<','<<kind<<','<<n<<','<<sample<<','<<it<<','<<std::setprecision(12)<<ns<<'\n';}

void benchmark(FftwApi& api,const Opt&o){
    if (o.n < 2 || !fftlab::pow2(o.n)) throw std::invalid_argument("vendor benchmark requires power-of-two N >= 2");
    if (o.samples < 5 || o.setup_samples < 1 || o.target_ms <= 0) throw std::invalid_argument("bad benchmark parameters");
    std::mt19937_64 rng(o.seed^o.n);std::uniform_real_distribution<double>d(-.5,.5);fftlab::Vector cx(o.n);for(auto&z:cx)z={d(rng),d(rng)};fftlab::RealVector rx(o.n);for(auto&v:rx)v=d(rng);
    fftlab::Radix2Plan ownc(o.n);fftlab::RealRadix2Plan ownr(o.n);auto ownbuf=cx;fftlab::HalfSpectrum hs(o.n/2+1);fftlab::Vector scratch(o.n/2);auto rbuf=rx;
    if (api.forget_wisdom) api.forget_wisdom();
    FftwComplexPair estc(api, o.n, FFTW_ESTIMATE);
    if (api.forget_wisdom) api.forget_wisdom();
    FftwRealPair estr(api, o.n, FFTW_ESTIMATE);
    estc.init(cx); estr.init(rx);
    if (api.forget_wisdom) api.forget_wisdom();
    FftwComplexPair measc(api, o.n, FFTW_MEASURE);
    if (api.forget_wisdom) api.forget_wisdom();
    FftwRealPair measr(api, o.n, FFTW_MEASURE);
    measc.init(cx); measr.init(rx);
    for(std::size_t w=0;w<o.warmups;++w){ownc.forward_inplace(ownbuf);ownc.inverse_inplace(ownbuf);ownr.forward(rbuf,hs,scratch);ownr.inverse(hs,rbuf,scratch);estc.pair();estr.pair();measc.pair();measr.pair();}
    auto pair_ownc=[&]{ownc.forward_inplace(ownbuf);ownc.inverse_inplace(ownbuf);}; auto it=calibrate(pair_ownc,o.target_ms);
    std::array<int,6> modes{0,1,2,3,4,5};
    if(o.raw)std::cout<<"phase,backend,planner,kind,N,sample,iterations,ns_per_operation\n";
    for(std::size_t s=0;s<o.samples;++s){std::shuffle(modes.begin(),modes.end(),rng);for(int m:modes){auto start=Clock::now();if(m==0){for(std::size_t i=0;i<it;++i)pair_ownc();}else if(m==1){for(std::size_t i=0;i<it;++i){ownr.forward(rbuf,hs,scratch);ownr.inverse(hs,rbuf,scratch);}}else if(m==2){for(std::size_t i=0;i<it;++i)estc.pair();}else if(m==3){for(std::size_t i=0;i<it;++i)estr.pair();}else if(m==4){for(std::size_t i=0;i<it;++i)measc.pair();}else{for(std::size_t i=0;i<it;++i)measr.pair();}auto ns=std::chrono::duration<double,std::nano>(Clock::now()-start).count()/double(2*it); if(m==0){sink=ownbuf[0].real();emit("execution","fftlab","native","complex",o.n,s,it,ns);}else if(m==1){sink=rbuf[0];emit("execution","fftlab","native","real",o.n,s,it,ns);}else if(m==2){sink=estc.data.p[0];emit("execution","fftw","estimate","complex",o.n,s,it,ns);}else if(m==3){sink=estr.real.p[0];emit("execution","fftw","estimate","real",o.n,s,it,ns);}else if(m==4){sink=measc.data.p[0];emit("execution","fftw","measure","complex",o.n,s,it,ns);}else{sink=measr.real.p[0];emit("execution","fftw","measure","real",o.n,s,it,ns);}}}
    // Setup samples. FFTW arrays are allocated outside the timed planner region, matching caller-buffer semantics.
    std::array<int,6> setup_modes{0,1,2,3,4,5};
    FftwBuffer ctmp(api,2*o.n),rtmp(api,o.n),stmp(api,2*(o.n/2+1));
    for(std::size_t s=0;s<o.setup_samples;++s){std::shuffle(setup_modes.begin(),setup_modes.end(),rng);for(int m:setup_modes){if((m>=2)&&api.forget_wisdom)api.forget_wisdom();auto start=Clock::now();if(m==0){fftlab::Radix2Plan p(o.n);sink=double(p.stored_twiddles());}else if(m==1){fftlab::RealRadix2Plan p(o.n);sink=double(p.spectrum_size());}else if(m==2||m==4){unsigned f=m==2?FFTW_ESTIMATE:FFTW_MEASURE;FftwPlanHandle pf(api,api.plan_dft(int(o.n),ctmp.p,ctmp.p,FFTW_FORWARD,f));FftwPlanHandle pb(api,api.plan_dft(int(o.n),ctmp.p,ctmp.p,FFTW_BACKWARD,f));sink=1;}else{unsigned f=m==3?FFTW_ESTIMATE:FFTW_MEASURE;FftwPlanHandle pf(api,api.plan_r2c(int(o.n),rtmp.p,stmp.p,f));FftwPlanHandle pb(api,api.plan_c2r(int(o.n),stmp.p,rtmp.p,f));sink=1;}auto ns=std::chrono::duration<double,std::nano>(Clock::now()-start).count();if(m==0)emit("setup","fftlab","native","complex",o.n,s,1,ns);else if(m==1)emit("setup","fftlab","native","real",o.n,s,1,ns);else if(m==2)emit("setup","fftw","estimate","complex",o.n,s,1,ns);else if(m==3)emit("setup","fftw","estimate","real",o.n,s,1,ns);else if(m==4)emit("setup","fftw","measure","complex",o.n,s,1,ns);else emit("setup","fftw","measure","real",o.n,s,1,ns);}}
}
}

int main(int argc,char**argv){try{auto o=parse(argc,argv);FftwApi api;if(o.info){std::cout<<"backend: fftw\navailable: "<<(api?"yes":"no")<<"\nlibrary: "<<api.library<<"\nversion: "<<api.version<<"\n";return api?0:77;}if(!api){std::cerr<<"FFTW runtime not found\n";return 77;}if(o.self){self_test(api);return 0;}if(o.raw){benchmark(api,o);return 0;}std::cerr<<"use --info, --self-test, or --raw-csv\n";return 2;}catch(const std::exception&e){std::cerr<<"error: "<<e.what()<<'\n';return 2;}}
