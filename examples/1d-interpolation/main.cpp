#include <iostream>
#include <fstream>
#include <pm.hpp>
#include <string>
#include <algorithm>

template<class callable>
void write(std::string fname, callable &F,int N_points)
{
    std::ofstream fd(fname+".dat");
    for(int i=0;i<N_points;++i)
    {
        double x=double(i)/N_points;
        fd << x << " " << F.interpolate({x}) << '\n';
    }
}

template<class callable>
void write_pts(std::string fname, callable &F)
{
    const int N_points = F.size();
    std::ofstream fd(fname+".dat");
    for(int i=0;i<N_points;++i)
    {
        double x=double(i)/N_points;
        fd << x << " " << F.interpolate({x}) << '\n';
    }
}

template <int k_max, class callable>
void test(int N_points, std::string name, callable F)
{
    constexpr int N = k_max*2+1;
    PM::grid<1,double, PM::Gaussian_filter, PM::Grid_filter<N> > exact(N);
    PM::grid<1,double, PM::Gaussian_filter, PM::NGP_filter > ngp(N);
    PM::grid<1,double, PM::Gaussian_filter, PM::CIC_filter > cic(N);
    PM::grid<1,double, PM::Gaussian_filter, PM::TSC_filter > tsc(N);
    PM::grid<1,double, PM::Gaussian_filter, PM::PCS_filter > pcs(N);
    PM::grid<1,double, PM::Gaussian_filter, PM::Gaussian_filter > gauss(N);
    PM::grid<1,double, PM::Gaussian_filter, PM::LowPass_filter<5,N> > low_pass(N);
    
    for(int i=0;i<N;++i)
    {
        exact.at({i})=ngp.at({i})=cic.at({i})=tsc.at({i})=pcs.at({i})=gauss.at({i})=
        low_pass.at({i})=
        F(double(i)/N);
    }
    
    write_pts("pts_"+name,exact);
    write("exact_"+name,exact,N_points);
    write("ngp_"+name,ngp,N_points);
    write("cic_"+name,cic,N_points);
    write("tsc_"+name,tsc,N_points);
    write("pcs_"+name,pcs,N_points);
    write("gauss_"+name,gauss,N_points);
    write("low_pass_"+name,low_pass,N_points);
}

template <int k_max, class callable>
void simple_test(int N_points, std::string name, callable F)
{
    constexpr int N = k_max*2+1;
    PM::grid<1,double, PM::Gaussian_filter, PM::LowPass_filter<k_max,N> > exact(N);
    for(int i=0;i<N;++i)
    {
        exact.at({i})=
        F(double(i)/N);
    }
    exact.interpolate({0});
}

const double pi = acos(-1);

int main()
{
    //simple_test<10>(10000,"triangle", [](double x) { return 3*std::min(x,1-x) ; });
    test<12>(10000,"triangle", [](double x) { return 3*std::min(x,1-x) ; });
    test<100>(10000,"sine_5_6", [](double x) { return sin(x*2*pi*6)+sin(x*2*pi*5); });
    return 0;
}
