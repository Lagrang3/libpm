#define BOOST_TEST_MODULE Field
#include <boost/mpi.hpp>
#include <boost/mpi/cartesian_communicator.hpp>
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

BOOST_AUTO_TEST_CASE(Size)
{
    mpi::cartesian_dimension dims[] = {{3, true}, {2, true}};
    BOOST_REQUIRE(fixture::world.size() == dims[0].size * dims[1].size);
    mpi::cartesian_communicator cart(fixture::world,
                                     mpi::cartesian_topology{dims});

    PM::Field<int> F(cart, 5);

    //  0   1
    //  ***|** 0
    //  ***|**
    //  ---+--
    //  ***|** 1
    //  ***|**
    //  ---+--
    //  ***|** 2

    size_t size{};
    switch (cart.rank())
    {
        case 0:
            size = 30;
            break;
        case 1:
            size = 20;
            break;
        case 2:
            size = 30;
            break;
        case 3:
            size = 20;
            break;
        case 4:
            size = 15;
            break;
        case 5:
            size = 10;
            break;
    }
    BOOST_CHECK_EQUAL(F.size(), size);
}
BOOST_AUTO_TEST_CASE(FFT)
{
    mpi::cartesian_dimension dims[] = {{3, true}, {2, true}};
    BOOST_REQUIRE(fixture::world.size() == dims[0].size * dims[1].size);
    mpi::cartesian_communicator cart(fixture::world,
                                     mpi::cartesian_topology{dims});

    const int N = 5;
    PM::Field<cd> F(cart, N);
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
BOOST_AUTO_TEST_CASE(transpose)
{
    mpi::cartesian_dimension dims[] = {{3, true}, {2, true}};
    BOOST_REQUIRE(fixture::world.size() == dims[0].size * dims[1].size);
    mpi::cartesian_communicator cart(fixture::world,
                                     mpi::cartesian_topology{dims});

    const int N = 5;
    PM::Field<std::array<int64_t, 3>> F(cart, N);

    //  0   1
    //  ***|** 0
    //  ***|**
    //  ---+--
    //  ***|** 1
    //  ***|**
    //  ---+--
    //  ***|** 2

    // setting the data
    const auto& offset = F.offset();
    for (int i = 0; i < F.extents()[0]; ++i)
        for (int j = 0; j < F.extents()[1]; ++j)
            for (int k = 0; k < F.extents()[2]; ++k)
                F(i, j, k) = {i + offset[0], j + offset[1], k + offset[2]};

    F.tranpose_yz();

    // check
    for (int i = 0; i < F.extents()[0]; ++i)
        for (int j = 0; j < F.extents()[1]; ++j)
            for (int k = 0; k < F.extents()[2]; ++k)
            {
                const auto& X = F(i, j, k);

                BOOST_CHECK_EQUAL(X[0], i + offset[0]);
                BOOST_CHECK_EQUAL(X[1], k + offset[2]);
                BOOST_CHECK_EQUAL(X[2], j + offset[1]);
            }

    F.tranpose_yz();
    for (int i = 0; i < F.extents()[0]; ++i)
        for (int j = 0; j < F.extents()[1]; ++j)
            for (int k = 0; k < F.extents()[2]; ++k)
            {
                const auto& X = F(i, j, k);

                BOOST_CHECK_EQUAL(X[0], i + offset[0]);
                BOOST_CHECK_EQUAL(X[1], j + offset[1]);
                BOOST_CHECK_EQUAL(X[2], k + offset[2]);
            }
    F.tranpose_xz();
    for (int i = 0; i < F.extents()[0]; ++i)
        for (int j = 0; j < F.extents()[1]; ++j)
            for (int k = 0; k < F.extents()[2]; ++k)
            {
                const auto& X = F(i, j, k);

                BOOST_CHECK_EQUAL(X[0], k + offset[2]);
                BOOST_CHECK_EQUAL(X[1], j + offset[1]);
                BOOST_CHECK_EQUAL(X[2], i + offset[0]);
            }
    F.tranpose_xz();
    for (int i = 0; i < F.extents()[0]; ++i)
        for (int j = 0; j < F.extents()[1]; ++j)
            for (int k = 0; k < F.extents()[2]; ++k)
            {
                const auto& X = F(i, j, k);

                BOOST_CHECK_EQUAL(X[0], i + offset[0]);
                BOOST_CHECK_EQUAL(X[1], j + offset[1]);
                BOOST_CHECK_EQUAL(X[2], k + offset[2]);
            }
}
