#define BOOST_TEST_MODULE Field Tranpose
#include <boost/mpi.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/environment.hpp>
#include <boost/test/data/monomorphic.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>
#include <pm/field.hpp>

#include <array>
#include <complex>
#include <pm/detail/utilities.hpp>

using cd = std::complex<double>;

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
    PM::Field<int> F(MPI_Comm(fixture::world), 5, {3, 2});

    //  0   1
    //  ***|** 0
    //  ***|**
    //  ---+--
    //  ***|** 1
    //  ***|**
    //  ---+--
    //  ***|** 2

    auto cart = F.get_communicator();

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
BOOST_AUTO_TEST_CASE(transpose)
{
    const int N = 5;
    PM::Field<std::array<int64_t, 3>> F(fixture::world, N, {3, 2});

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
