#include <algorithm>
#include <charconv>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <numbers>
#include <random>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace fftlab {
using Complex = std::complex<double>;
using Vector = std::vector<Complex>;
constexpr double pi = std::numbers::pi_v<double>;

bool pow2(std::size_t n) { return n && !(n & (n - 1)); }
std::size_t next_pow2(std::size_t n) {
    if (n <= 1) return 1;
    const auto high = std::size_t{1} << (std::numeric_limits<std::size_t>::digits - 1);
    if (n > high) throw std::overflow_error("next power of two overflows size_t");
    --n;
    for (std::size_t s = 1; s < std::numeric_limits<std::size_t>::digits; s <<= 1) n |= n >> s;
    return n + 1;
}
Complex root(double a) { return {std::cos(a), std::sin(a)}; }

Vector dft(const Vector& x, bool inv = false) {
    const auto n = x.size(); Vector y(n);
    if (!n) return y;
    const double sign = inv ? 1.0 : -1.0, scale = inv ? 1.0 / double(n) : 1.0;
    for (std::size_t k = 0; k < n; ++k) {
        Complex s{};
        for (std::size_t t = 0; t < n; ++t)
            s += x[t] * root(sign * 2.0 * pi * double(k) * double(t) / double(n));
        y[k] = s * scale;
    }
    return y;
}

void radix2_inplace(Vector& a, bool inv = false) {
    const auto n = a.size();
    if (!n) return;
    if (!pow2(n)) throw std::invalid_argument("radix2 requires power-of-two N");
    for (std::size_t i = 1, j = 0; i < n; ++i) {
        std::size_t bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) std::swap(a[i], a[j]);
    }
    const double sign = inv ? 1.0 : -1.0;
    for (std::size_t len = 2; len <= n;) {
        const Complex step = root(sign * 2.0 * pi / double(len));
        for (std::size_t i = 0; i < n; i += len) {
            Complex w{1, 0};
            for (std::size_t j = 0; j < len / 2; ++j) {
                const Complex u = a[i+j], v = a[i+j+len/2] * w;
                a[i+j] = u + v; a[i+j+len/2] = u - v; w *= step;
            }
        }
        if (len == n) break;
        len <<= 1;
    }
    if (inv) for (auto& z : a) z /= double(n);
}
Vector radix2(const Vector& x, bool inv = false) { Vector y = x; radix2_inplace(y, inv); return y; }

void recursive_core(Vector& a, bool inv) {
    const auto n = a.size(); if (n <= 1) return;
    Vector e(n/2), o(n/2);
    for (std::size_t i = 0; i < n/2; ++i) { e[i] = a[2*i]; o[i] = a[2*i+1]; }
    recursive_core(e, inv); recursive_core(o, inv);
    Complex w{1,0}, step = root((inv ? 1.0 : -1.0) * 2.0 * pi / double(n));
    for (std::size_t k = 0; k < n/2; ++k) {
        const auto t = w * o[k]; a[k] = e[k] + t; a[k+n/2] = e[k] - t; w *= step;
    }
}
Vector recursive(const Vector& x, bool inv = false) {
    if (!x.empty() && !pow2(x.size())) throw std::invalid_argument("recursive radix2 requires power-of-two N");
    Vector y = x; recursive_core(y, inv); if (inv && !y.empty()) for (auto& z : y) z /= double(y.size()); return y;
}

std::size_t factor(std::size_t n) {
    if (!(n % 2)) return 2;
    for (std::size_t p = 3; p <= n/p; p += 2) if (!(n % p)) return p;
    return n;
}
Vector mixed_core(const Vector& x, bool inv) {
    const auto n = x.size(); if (n <= 1) return x;
    const auto r = factor(n);
    if (r == n) {
        Vector y(n); const double sign = inv ? 1.0 : -1.0;
        for (std::size_t k=0;k<n;++k) for (std::size_t t=0;t<n;++t)
            y[k] += x[t] * root(sign * 2.0*pi*double(k)*double(t)/double(n));
        return y;
    }
    const auto m = n/r; std::vector<Vector> sub(r, Vector(m));
    for (std::size_t q=0;q<r;++q) {
        for (std::size_t j=0;j<m;++j) sub[q][j]=x[r*j+q];
        sub[q]=mixed_core(sub[q],inv);
    }
    Vector y(n); const double sign = inv ? 1.0 : -1.0;
    for (std::size_t k0=0;k0<m;++k0) for (std::size_t k1=0;k1<r;++k1) {
        const auto k=k0+m*k1; Complex w{1,0}, step=root(sign*2.0*pi*double(k)/double(n));
        for (std::size_t q=0;q<r;++q) { y[k]+=sub[q][k0]*w; w*=step; }
    }
    return y;
}
Vector mixed(const Vector& x, bool inv=false) {
    Vector y=mixed_core(x,inv); if(inv&&!y.empty()) for(auto& z:y) z/=double(y.size()); return y;
}

Vector bluestein(const Vector& x, bool inv=false) {
    const auto n=x.size(); if(n<=1) return x;
    if(n>(std::numeric_limits<std::size_t>::max()/2)+1) throw std::length_error("Bluestein workspace overflow");
    const auto m=next_pow2(2*n-1); Vector a(m),b(m); const double sign=inv?1.0:-1.0;
    for(std::size_t k=0;k<n;++k){
        const long double kd=k, period=2.0L*n, phase=std::fmod(kd*kd,period)/n;
        const Complex c=root(sign*pi*double(phase)); a[k]=x[k]*c; b[k]=std::conj(c); if(k) b[m-k]=std::conj(c);
    }
    radix2_inplace(a); radix2_inplace(b); for(std::size_t i=0;i<m;++i)a[i]*=b[i]; radix2_inplace(a,true);
    Vector y(n);
    for(std::size_t k=0;k<n;++k){
        const long double kd=k, period=2.0L*n, phase=std::fmod(kd*kd,period)/n;
        y[k]=a[k]*root(sign*pi*double(phase)); if(inv)y[k]/=double(n);
    }
    return y;
}

bool smooth235(std::size_t n){ if(!n)return false; for(auto p:{2u,3u,5u})while(!(n%p))n/=p; return n==1; }
enum class Algo { Auto,Dft,Radix2,Recursive,Mixed,Bluestein };
std::string_view name(Algo a){ switch(a){case Algo::Auto:return"auto";case Algo::Dft:return"dft";case Algo::Radix2:return"radix2-iterative";case Algo::Recursive:return"radix2-recursive";case Algo::Mixed:return"mixed-radix";case Algo::Bluestein:return"bluestein";} return"?"; }
Algo parse_algo(std::string_view s){ if(s=="auto")return Algo::Auto;if(s=="dft")return Algo::Dft;if(s=="radix2"||s=="radix2-iterative")return Algo::Radix2;if(s=="recursive"||s=="radix2-recursive")return Algo::Recursive;if(s=="mixed"||s=="mixed-radix")return Algo::Mixed;if(s=="bluestein"||s=="chirpz")return Algo::Bluestein;throw std::invalid_argument("unknown algorithm"); }
bool supports(Algo a,std::size_t n){return (a!=Algo::Radix2&&a!=Algo::Recursive)||!n||pow2(n);}
Vector transform(const Vector& x,Algo a,bool inv=false){
    if(a==Algo::Auto){ if(x.size()<=1)return x; if(pow2(x.size()))return radix2(x,inv); if(smooth235(x.size()))return mixed(x,inv); return bluestein(x,inv); }
    if(a==Algo::Dft)return dft(x,inv); if(a==Algo::Radix2)return radix2(x,inv); if(a==Algo::Recursive)return recursive(x,inv); if(a==Algo::Mixed)return mixed(x,inv); return bluestein(x,inv);
}
Vector signal(std::size_t n){ std::mt19937_64 g(0xF17F17ULL^n);std::uniform_real_distribution<double>d(-.5,.5);Vector x(n);for(std::size_t i=0;i<n;++i){double t=n?double(i)/n:0;x[i]={.9*std::sin(2*pi*3*t)+.35*std::cos(2*pi*11*t)+.05*d(g),.4*std::sin(2*pi*5*t)+.05*d(g)};}return x; }

using LComplex=std::complex<long double>; using LVector=std::vector<LComplex>;
LVector oracle(const Vector& x){ const auto n=x.size();LVector y(n);constexpr long double p=std::numbers::pi_v<long double>;for(std::size_t k=0;k<n;++k)for(std::size_t t=0;t<n;++t)y[k]+=LComplex{x[t].real(),x[t].imag()}*LComplex{std::cos(-2*p*k*t/n),std::sin(-2*p*k*t/n)};return y; }
struct Err{double max{},rms{},rel{};};
Err error(const Vector&a,const LVector&b){long double mx=0,se=0,sr=0;for(std::size_t i=0;i<a.size();++i){LComplex z{a[i].real(),a[i].imag()};auto e=std::abs(z-b[i]);mx=std::max(mx,e);se+=e*e;auto r=std::abs(b[i]);sr+=r*r;}long double n=a.empty()?1:a.size();return{double(mx),double(std::sqrt(se/n)),double(sr?std::sqrt(se/sr):std::sqrt(se))};}
double roundtrip(const Vector&x,Algo a){auto y=transform(transform(x,a),a,true);double e=0;for(std::size_t i=0;i<x.size();++i)e=std::max(e,std::abs(y[i]-x[i]));return e;}

volatile double sink=0;
struct Stats{std::size_t n{},samples{},iters{};Algo algo{};double min{},p05{},median{},mean{},p95{},max{},sd{};};
double pct(const std::vector<double>&v,double p){double pos=p*(v.size()-1);auto lo=std::size_t(std::floor(pos)),hi=std::size_t(std::ceil(pos));return lo==hi?v[lo]:v[lo]*(hi-pos)+v[hi]*(pos-lo);}
std::size_t calibrate(const Vector&x,Algo a,double ms){using C=std::chrono::steady_clock;std::size_t it=1;while(it<(1u<<20)){auto s=C::now();double c=0;for(std::size_t i=0;i<it;++i){auto y=transform(x,a);if(!y.empty())c+=y[i%y.size()].real();}sink=c;double e=std::chrono::duration<double,std::milli>(C::now()-s).count();if(e>=ms*.5)break;it=std::min<std::size_t>(1u<<20,std::max(it+1,std::size_t(it*std::clamp(ms/std::max(e,1e-9),2.0,16.0))));}return it;}
Stats bench(const Vector&x,Algo a,std::size_t samples,std::size_t warm,double target){if(!supports(a,x.size()))throw std::invalid_argument("unsupported size");if(samples<3||target<=0)throw std::invalid_argument("bad benchmark parameters");for(std::size_t i=0;i<warm;++i){auto y=transform(x,a);if(!y.empty())sink=y[0].real();}auto it=calibrate(x,a,target);std::vector<double>v;using C=std::chrono::steady_clock;for(std::size_t s=0;s<samples;++s){auto st=C::now();double c=0;for(std::size_t i=0;i<it;++i){auto y=transform(x,a);if(!y.empty())c+=y[(i+s)%y.size()].real();}sink=c;v.push_back(std::chrono::duration<double,std::nano>(C::now()-st).count()/it);}std::sort(v.begin(),v.end());double mean=std::accumulate(v.begin(),v.end(),0.0)/v.size(),var=0;for(double z:v)var+=(z-mean)*(z-mean);var/=v.size()-1;return{x.size(),samples,it,a,v.front(),pct(v,.05),pct(v,.5),mean,pct(v,.95),v.back(),std::sqrt(var)};}

void tests(){std::size_t checks=0;auto req=[&](bool c){++checks;if(!c)throw std::runtime_error("self-test failure");};req(next_pow2(3)==4);for(std::size_t n:{2u,3u,5u,6u,7u,8u,10u,12u,16u,25u,64u}){auto x=signal(n),ref=dft(x);for(auto a:{Algo::Auto,Algo::Mixed,Algo::Bluestein}){auto y=transform(x,a);for(std::size_t i=0;i<n;++i)req(std::abs(y[i]-ref[i])<2e-9*(1+std::abs(ref[i])));req(roundtrip(x,a)<2e-9);}if(pow2(n))for(auto a:{Algo::Radix2,Algo::Recursive}){auto y=transform(x,a);for(std::size_t i=0;i<n;++i)req(std::abs(y[i]-ref[i])<2e-9*(1+std::abs(ref[i])));req(roundtrip(x,a)<2e-9);}}bool threw=false;try{(void)radix2(signal(6));}catch(const std::invalid_argument&){threw=true;}req(threw);std::cout<<"PASS: "<<checks<<" checks\n";}

struct Opt{bool help=false,self=false,verify=false,one=false,suite=false,csv=false;Algo algo=Algo::Auto;std::size_t size=1024,samples=31,warm=5;double target=5;std::vector<std::size_t>sizes{8,16,32,64,128,256,512,1000,1009,1024,2048,4096};};
std::size_t number(std::string_view s){std::size_t n{};auto[p,e]=std::from_chars(s.data(),s.data()+s.size(),n);if(e!=std::errc{}||p!=s.data()+s.size()||!n)throw std::invalid_argument("invalid integer");return n;}
std::vector<std::size_t> list(std::string_view s){std::vector<std::size_t>v;while(!s.empty()){auto p=s.find(',');v.push_back(number(s.substr(0,p)));if(p==s.npos)break;s.remove_prefix(p+1);}return v;}
Opt options(int ac,char**av){Opt o;auto val=[&](int&i){if(++i>=ac)throw std::invalid_argument("missing option value");return std::string_view(av[i]);};for(int i=1;i<ac;++i){std::string_view a=av[i];if(a=="-h"||a=="--help")o.help=true;else if(a=="--self-test")o.self=true;else if(a=="--verify")o.verify=true;else if(a=="--benchmark")o.one=true;else if(a=="--benchmark-suite")o.suite=true;else if(a=="--csv")o.csv=true;else if(a=="--algorithm")o.algo=parse_algo(val(i));else if(a=="--size")o.size=number(val(i));else if(a=="--samples")o.samples=number(val(i));else if(a=="--warmups")o.warm=number(val(i));else if(a=="--target-ms")o.target=std::stod(std::string(val(i)));else if(a=="--sizes")o.sizes=list(val(i));else throw std::invalid_argument("unknown option");}return o;}
void help(){std::cout<<"fft - dependency-free C++23 Fourier transform laboratory\n\n--self-test\n--verify --algorithm ALGO --size N\n--benchmark --algorithm ALGO --size N [--samples 31] [--target-ms 5] [--csv]\n--benchmark-suite [--sizes 8,16,...] [--csv]\n\nAlgorithms: auto, dft, radix2-iterative, radix2-recursive, mixed-radix, bluestein\n";}
void csv_head(){std::cout<<"algorithm,N,samples,iterations_per_sample,min_ns,p05_ns,median_ns,mean_ns,p95_ns,max_ns,stddev_ns\n";}
void print(const Stats&s,bool csv){if(csv)std::cout<<name(s.algo)<<','<<s.n<<','<<s.samples<<','<<s.iters<<','<<std::setprecision(12)<<s.min<<','<<s.p05<<','<<s.median<<','<<s.mean<<','<<s.p95<<','<<s.max<<','<<s.sd<<'\n';else std::cout<<std::left<<std::setw(18)<<name(s.algo)<<" N="<<std::setw(6)<<s.n<<" median="<<std::fixed<<std::setprecision(1)<<s.median<<" ns p05="<<s.p05<<" p95="<<s.p95<<" sd="<<s.sd<<'\n';}
std::vector<Algo> suite(std::size_t n){std::vector<Algo>v{Algo::Auto,Algo::Mixed,Algo::Bluestein};if(pow2(n)){v.push_back(Algo::Radix2);v.push_back(Algo::Recursive);}if(n<=2048)v.push_back(Algo::Dft);return v;}
}

int main(int argc,char**argv){using namespace fftlab;try{auto o=options(argc,argv);if(argc==1||o.help){help();return 0;}if(o.self){tests();return 0;}if(o.verify){if(o.size>4096)throw std::invalid_argument("verification capped at N=4096");auto x=signal(o.size);auto e=error(transform(x,o.algo),oracle(x));std::cout<<std::setprecision(12)<<"algorithm: "<<name(o.algo)<<"\nN: "<<o.size<<"\nmax_abs_error: "<<e.max<<"\nrms_error: "<<e.rms<<"\nrelative_l2_error: "<<e.rel<<"\nround_trip_max_abs_error: "<<roundtrip(x,o.algo)<<'\n';return 0;}if(o.one){if(o.csv)csv_head();print(bench(signal(o.size),o.algo,o.samples,o.warm,o.target),o.csv);return 0;}if(o.suite){if(o.csv)csv_head();for(auto n:o.sizes)for(auto a:suite(n))print(bench(signal(n),a,o.samples,o.warm,o.target),o.csv);return 0;}throw std::invalid_argument("no action requested");}catch(const std::exception&e){std::cerr<<"error: "<<e.what()<<'\n';return 2;}}
