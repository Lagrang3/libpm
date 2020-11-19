/*
    Unit Test:
    ----------

    This test module checks the interpolation
    routine of the Field class.

    Should this test fail the possible causes are:

    1. the grid size N is too small for the interpolation
    to reproduce the input function, for the given
    tolerance. If the input function
    has no moments higher than N/2 and the
    interpolating kernel is the Shannon-Whittaker,
    here PM::Grid_filter, then interpolation should
    be exact.

    2. PM::Field::interpolate is buggy.
*/

#define BOOST_TEST_MODULE 1D Interpolation
#include <boost/mpi.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/environment.hpp>
#include <boost/mpl/list.hpp>
#include <boost/test/unit_test.hpp>
#include <cmath>
#include <limits>
#include <pm/field.hpp>
#include <pm/filters.hpp>
#include <vector>

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

template <class filter_t, class callable>
double interpolate_n_check(const int N_grid, const int N_eval, callable myfun)
{
    PM::Field<double, filter_t> F(fixture::world, N_grid, {2, 1});

    // const auto& offset = F.offset();
    for (int i = 0; i < F.extents()[0]; ++i)
        for (int j = 0; j < F.extents()[1]; ++j)
            for (int k = 0; k < F.extents()[2]; ++k)
                F(i, j, k) = myfun(k * 1.0 / N_grid);

    F.update_ghosts();

    double max_diff = 0;
    // for (int i = 0; i < F.extents()[0]; ++i)
    //    for (int j = 0; j < F.extents()[1]; ++j)
    {
        int i = 0, j = 0;
        for (int k = 0; k < N_eval; ++k)
        {
            double x = i, y = j, z = k * (N_grid * 1.0 / N_eval);
            // double diff = std::abs(F.interpolate(x,y,z) - myfun(z/N_grid));
            double fi = F.interpolate(x, y, z), fo = myfun(z / N_grid);
            double diff = std::abs(fi - fo);
            // if(fixture::world.rank()==0)std::cout << z << " " <<  fi << " "
            // << fo << "\n";
            max_diff = std::max(max_diff, diff);
        }
    }
    return max_diff;
}

template <class filter_t, class callable>
double check_homogenuity(const int N_grid, const int N_eval, callable myfun)
{
    PM::Field<double, filter_t> F(fixture::world, N_grid, {2, 1});

    // const auto& offset = F.offset();
    for (int i = 0; i < F.extents()[0]; ++i)
        for (int j = 0; j < F.extents()[1]; ++j)
            for (int k = 0; k < F.extents()[2]; ++k)
                F(i, j, k) = myfun(k * 1.0 / N_grid);

    F.update_ghosts();

    std::vector<double> values(N_eval);
    {
        int i = 0, j = 0;
        for (int k = 0; k < N_eval; ++k)
        {
            double x = i, y = j, z = k * (N_grid * 1.0 / N_eval);
            values[k] = F.interpolate(x, y, z);
        }
    }

    double max_diff = 0;
    for (int i = 0; i < F.extents()[0]; ++i)
        for (int j = 0; j < F.extents()[1]; ++j)
            for (int k = 0; k < N_eval; ++k)
            {
                double x = i, y = j, z = k * (N_grid * 1.0 / N_eval);
                double fi = F.interpolate(x, y, z), fo = values[k];
                double diff = std::abs(fi - fo);
                max_diff = std::max(max_diff, diff);
            }
    return max_diff;
}

/*
    This tests the interpolation accuracy of several kernels, for a couple of
    1-dimensional input functions. The `Field` is still 3-dimensional, so the
    field is constant along x and y axis.
*/
BOOST_AUTO_TEST_CASE(oned_accuracy)
{
    auto Triangle_shape = [](double z) { return 3 * std::min(z, 1 - z); };
    auto Sine_shape = [](double z) {
        static const double pi = acos(-1.0);
        return sin(2 * pi * z) + 3 * sin(2 * pi * z * 2);
    };

    const int N_grid = 20;
    const int N_eval = 200;

    //  0
    //  ***** 0
    //  *****
    //  *****
    //  -----
    //  ***** 1
    //  *****
    double d;

    d = interpolate_n_check<PM::filters::CIC>(N_grid, N_eval, Triangle_shape);
    BOOST_CHECK_SMALL(d, 0.1);
    d = interpolate_n_check<PM::filters::CIC>(N_grid, N_eval, Sine_shape);
    BOOST_CHECK_SMALL(d, 0.2);

    d = interpolate_n_check<PM::filters::NGP>(N_grid, N_eval, Triangle_shape);
    BOOST_CHECK_SMALL(d, 0.2);
    d = interpolate_n_check<PM::filters::NGP>(N_grid, N_eval, Sine_shape);
    BOOST_CHECK_SMALL(d, 1.0);

    d = interpolate_n_check<PM::filters::TSC>(N_grid, N_eval, Triangle_shape);
    BOOST_CHECK_SMALL(d, 0.2);
    d = interpolate_n_check<PM::filters::TSC>(N_grid, N_eval, Sine_shape);
    BOOST_CHECK_SMALL(d, 0.2);

    d = interpolate_n_check<PM::filters::PCS>(N_grid, N_eval, Triangle_shape);
    BOOST_CHECK_SMALL(d, 0.2);
    d = interpolate_n_check<PM::filters::PCS>(N_grid, N_eval, Sine_shape);
    BOOST_CHECK_SMALL(d, 0.3);

    d = interpolate_n_check<PM::filters::Gaussian>(N_grid, N_eval,
                                                   Triangle_shape);
    BOOST_CHECK_SMALL(d, 0.2);
    d = interpolate_n_check<PM::filters::Gaussian>(N_grid, N_eval, Sine_shape);
    BOOST_CHECK_SMALL(d, 1.0);
}
/*
    The `Field` is still 3-dimensional and it is being tested with a
    1-dimensional input function that runs along the z axis, the rest remain
    constant. This tests that the interpolation is indeed homogeneous along x
    and y.
*/
BOOST_AUTO_TEST_CASE(homegenuity)
{
    auto Triangle_shape = [](double z) { return 3 * std::min(z, 1 - z); };
    auto Sine_shape = [](double z) {
        static const double pi = acos(-1.0);
        return sin(2 * pi * z) + 3 * sin(2 * pi * z * 2);
    };

    const int N_grid = 20;
    const int N_eval = 200;

    //  0
    //  ***** 0
    //  *****
    //  *****
    //  -----
    //  ***** 1
    //  *****
    double d;
    const double epsilon = std::numeric_limits<double>::epsilon();

    d = check_homogenuity<PM::filters::CIC>(N_grid, N_eval, Triangle_shape);
    BOOST_CHECK_SMALL(d, epsilon);
    d = check_homogenuity<PM::filters::CIC>(N_grid, N_eval, Sine_shape);
    BOOST_CHECK_SMALL(d, epsilon);

    d = check_homogenuity<PM::filters::NGP>(N_grid, N_eval, Triangle_shape);
    BOOST_CHECK_SMALL(d, epsilon);
    d = check_homogenuity<PM::filters::NGP>(N_grid, N_eval, Sine_shape);
    BOOST_CHECK_SMALL(d, epsilon);

    d = check_homogenuity<PM::filters::TSC>(N_grid, N_eval, Triangle_shape);
    BOOST_CHECK_SMALL(d, epsilon);
    d = check_homogenuity<PM::filters::TSC>(N_grid, N_eval, Sine_shape);
    BOOST_CHECK_SMALL(d, epsilon);

    d = check_homogenuity<PM::filters::PCS>(N_grid, N_eval, Triangle_shape);
    BOOST_CHECK_SMALL(d, epsilon);
    d = check_homogenuity<PM::filters::PCS>(N_grid, N_eval, Sine_shape);
    BOOST_CHECK_SMALL(d, epsilon);

    d = check_homogenuity<PM::filters::Gaussian>(N_grid, N_eval,
                                                 Triangle_shape);
    BOOST_CHECK_SMALL(d, epsilon);
    d = check_homogenuity<PM::filters::Gaussian>(N_grid, N_eval, Sine_shape);
    BOOST_CHECK_SMALL(d, epsilon);
}
