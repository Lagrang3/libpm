#define BOOST_TEST_MODULE Field FFT
#include <boost/mpi.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/environment.hpp>
#include <boost/test/data/monomorphic.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>
#include <pm/field.hpp>

#include <array>
#include <complex>
#include <fftw3.h>
#include <pm/detail/utilities.hpp>

using cd = std::complex<double>;

/*
    From Stroustrup's The C++ Programming Language 4ed
    section 13.3.1
*/
template <class F>
struct Final_action
{
    F clean;
    Final_action(F f) : clean{f} {}
    ~Final_action() { clean(); }
};
template <class F>
Final_action<F> finally(F f)
{
    return Final_action<F>(f);
}

namespace ut = boost::unit_test;
namespace mpi = boost::mpi;

struct fixture
{
    fixture() {}

    ~fixture() {}

    static mpi::environment env;
    static mpi::communicator world;
};

mpi::environment fixture::env;
mpi::communicator fixture::world;

BOOST_TEST_GLOBAL_FIXTURE(fixture);

BOOST_AUTO_TEST_CASE(FFT)
{
    const int N = 5;
    PM::Field<cd> F(fixture::world, N, {3, 2});
    std::array<int64_t, 3> N_strides{N * N, N, 1};

    //  0   1
    //  ***|** 0
    //  ***|**
    //  ---+--
    //  ***|** 1
    //  ***|**
    //  ---+--
    //  ***|** 2

    // global function
    auto myfun = [](int x, int y, int z) -> double {
        double res = x * sqrt(1.0 * y) * sin(z * 1.0 / N);
        return res;
    };

    // fftw initialization
    fftw_complex* in = new fftw_complex[N * N * N];
    fftw_plan p =
        fftw_plan_dft_3d(N, N, N, in, in, FFTW_FORWARD, FFTW_ESTIMATE);
    auto destroy_plan = finally([&]() {
        delete[] in;
        fftw_destroy_plan(p);
    });

    // setting the data
    const auto& offset = F.offset();
    for (int i = 0; i < F.extents()[0]; ++i)
        for (int j = 0; j < F.extents()[1]; ++j)
            for (int k = 0; k < F.extents()[2]; ++k)
                F(i, j, k) = myfun(i + offset[0], j + offset[1], k + offset[2]);

    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            for (int k = 0; k < N; ++k)
            {
                in[PM::utilities::index(i, j, k, N_strides)][0] =
                    myfun(i, j, k);
                in[PM::utilities::index(i, j, k, N_strides)][1] = 0.0;
            }

    // execute
    F.fft();
    fftw_execute(p);

    // compare
    double diff = 0;
    for (int i = 0; i < F.extents()[0]; ++i)
        for (int j = 0; j < F.extents()[1]; ++j)
            for (int k = 0; k < F.extents()[2]; ++k)
            {
                fftw_complex& v = in[PM::utilities::index(
                    i + offset[0], j + offset[1], k + offset[2], N_strides)];
                diff += std::abs(F(i, j, k) - cd{v[0], v[1]});
            }
    BOOST_CHECK_SMALL(diff, 1e-12);
}
