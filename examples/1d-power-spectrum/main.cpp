#include <fstream>
#include <iostream>
#include <pm.hpp>
#include <random>
#include <string>
#include <vector>

void write(std::vector<double>& modes, std::string fname)
{
    std::ofstream fd("pw-random_" + fname + ".dat");
    for (auto x : modes)
        fd << x << " ";
    fd << '\n';
}

auto random_positions(int Np /* # particles */)
{
    std::vector<double> pos;
    std::mt19937 rng;
    std::normal_distribution<double> dist(10, 0.005);
    for (int i = 0; i < Np; ++i)
        pos.push_back(dist(rng));
    return pos;
}

template <int N, class filter_t>
void power_spectrum(const std::vector<double>& pos, std::string fname)
{
    PM::grid<1, double, filter_t, PM::Grid_filter<N>> mygrid(N);
    mygrid.sample_density(pos);
    mygrid.fft();
    auto modes = mygrid.get_modes();
    double fact = 1.0 / pos.size();
    fact *= fact;
    for (auto& x : modes)
        x *= fact;
    write(modes, fname);
}

int main()
{
    auto pos = random_positions(10000);
    power_spectrum<100, PM::Grid_filter<100>>(pos, "exact");
    power_spectrum<100, PM::NGP_filter>(pos, "ngp");
    power_spectrum<100, PM::CIC_filter>(pos, "cic");
    power_spectrum<100, PM::TSC_filter>(pos, "tsc");
    power_spectrum<100, PM::PCS_filter>(pos, "pcs");
    power_spectrum<100, PM::Gaussian_filter>(pos, "gauss");
    power_spectrum<100, PM::LowPass_filter<10, 100>>(pos, "low_pass");
    return 0;
}
