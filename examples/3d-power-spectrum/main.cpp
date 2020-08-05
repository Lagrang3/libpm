#include <chrono>
#include <fstream>
#include <gadget.hpp>
#include <iostream>
#include <pm.hpp>
#include <random>
#include <string>
#include <vector>

template <class T>
void write(std::vector<T>& modes, std::string fname)
{
    std::ofstream fd("pw-gevolution_" + fname + ".dat");
    for (auto x : modes)
        fd << x << " ";
    fd << '\n';
}

auto get_positions()
{
    gadget::isnapshot<1> snap("newton_snap006_cdm");
    int npart = 0;
    for (int i = 0; i < gadget::PTYPES; ++i)
        npart += snap.npart(i);

    std::vector<float> buff(3 * npart);

    snap.scan_block("POS ", buff.begin());
    float ilength = 1 / snap.header.get_BoxSize();
    for (auto& x : buff)
        x *= ilength;
    // return std::vector<double>(buff.begin(),buff.end());
    return buff;
}
struct timeit
{
    std::string message;
    std::chrono::steady_clock::time_point start;

    timeit(std::string mes)
        : message(mes), start(std::chrono::steady_clock::now())
    {
    }
    ~timeit()
    {
        auto dt = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::steady_clock::now() - start);
        std::cerr << "Timeit: " << message << " --> " << dt.count() << " ms\n";
    }
};

template <int N, class filter_t>
void power_spectrum(const std::vector<float>& pos, std::string fname)
{
    std::cout << "Test Case " << fname << "\n";
    PM::grid<3, float, filter_t, PM::Grid_filter<N>> mygrid(N);
    {
        timeit T("sampling density");
        mygrid.sample_density(pos);
    }
    {
        timeit T("performing fft");
        mygrid.fft();
    }
    {
        timeit T("performing mode correction");
        mygrid.sample_correction();
    }
    std::vector<float> modes;
    {
        timeit T("evaluating the modes");
        modes = mygrid.get_modes();
    }
    float fact = 3.0 / pos.size(), two_pi = 2 * acos(-1.0);
    fact *= fact / two_pi / two_pi / two_pi;

    for (auto& x : modes)
        x *= fact;
    write(modes, fname);
}

int main()
{
    std::vector<float> pos;
    {
        timeit T("reading snapshot");
        pos = get_positions();
    }
    power_spectrum<256, PM::NGP_filter>(pos, "ngp");
    power_spectrum<256, PM::CIC_filter>(pos, "cic");
    power_spectrum<256, PM::TSC_filter>(pos, "tsc");
    power_spectrum<256, PM::PCS_filter>(pos, "pcs");
    power_spectrum<256, PM::Gaussian_filter>(pos, "gauss");
    return 0;
}
