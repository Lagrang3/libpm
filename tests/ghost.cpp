#define BOOST_TEST_MODULE Ghost
#include <boost/mpi.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/environment.hpp>
#include <boost/test/data/monomorphic.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>
#include <pm/field.hpp>

#include <array>
#include <complex>

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

BOOST_AUTO_TEST_CASE(ghost_indexing)
{
    using PM::Field;
    using PM::utilities::modulo;

    const int N = 5;
    Field<std::array<int64_t, 3>> F(fixture::world, N, {3, 2});

    //  0   1
    //  ***|** 0
    //  ***|**
    //  ---+--
    //  ***|** 1
    //  ***|**
    //  ---+--
    //  ***|** 2

    const int nghost = F.ghost_thick();

    BOOST_TEST_MESSAGE("x,y,z in-range");
    for (int i = 0; i < F.extents()[0]; ++i)
        for (int j = 0; j < F.extents()[1]; ++j)
            for (int k = 0; k < F.extents()[2]; ++k)
                BOOST_CHECK_NO_THROW(F.at(i, j, k));

    BOOST_TEST_MESSAGE("z off-range");
    for (int i = 0; i < F.extents()[0]; ++i)
        for (int j = 0; j < F.extents()[1]; ++j)
            for (int k = -nghost; k < F.extents()[2] + nghost; ++k)
                BOOST_CHECK_NO_THROW(F.at(i, j, k));

    BOOST_TEST_MESSAGE("y off-range");
    for (int i = 0; i < F.extents()[0]; ++i)
        for (int j = -nghost; j < F.extents()[1] + nghost; ++j)
            for (int k = 0; k < F.extents()[2]; ++k)
                BOOST_CHECK_NO_THROW(F.at(i, j, k));

    BOOST_TEST_MESSAGE("x off-range");
    for (int i = -nghost; i < F.extents()[0] + nghost; ++i)
        for (int j = 0; j < F.extents()[1]; ++j)
            for (int k = 0; k < F.extents()[2]; ++k)
                BOOST_CHECK_NO_THROW(F.at(i, j, k));

    BOOST_TEST_MESSAGE("y,z off-range");
    for (int i = 0; i < F.extents()[0]; ++i)
        for (int j = 0 - nghost; j < F.extents()[1] + nghost; ++j)
            for (int k = -nghost; k < F.extents()[2] + nghost; ++k)
                BOOST_CHECK_NO_THROW(F.at(i, j, k));

    BOOST_TEST_MESSAGE("x,z off-range");
    for (int i = 0 - nghost; i < F.extents()[0] + nghost; ++i)
        for (int j = 0; j < F.extents()[1]; ++j)
            for (int k = -nghost; k < F.extents()[2] + nghost; ++k)
                BOOST_CHECK_NO_THROW(F.at(i, j, k));

    BOOST_TEST_MESSAGE("x,y off-range");
    for (int i = 0 - nghost; i < F.extents()[0] + nghost; ++i)
        for (int j = -nghost; j < F.extents()[1] + nghost; ++j)
            for (int k = 0; k < F.extents()[2]; ++k)
                BOOST_CHECK_NO_THROW(F.at(i, j, k));

    BOOST_TEST_MESSAGE("x,y,z off-range");
    for (int i = -nghost; i < F.extents()[0] + nghost; ++i)
        for (int j = -nghost; j < F.extents()[1] + nghost; ++j)
            for (int k = -nghost; k < F.extents()[2] + nghost; ++k)
                BOOST_CHECK_NO_THROW(F.at(i, j, k));
}

BOOST_AUTO_TEST_CASE(ghost_update)
{
    using PM::Field;
    using PM::utilities::modulo;

    const int N = 5;
    Field<std::array<int64_t, 3>> F(fixture::world, N, {3, 2});

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

    F.update_ghosts();

    const int nghost = F.ghost_thick();
    for (int i = -nghost; i < F.extents()[0] + nghost; ++i)
        for (int j = -nghost; j < F.extents()[1] + nghost; ++j)
            for (int k = -nghost; k < F.extents()[2] + nghost; ++k)
            {
                const auto& X = F(i, j, k);

                BOOST_CHECK_EQUAL(X[0], modulo(i + offset[0], N));
                BOOST_CHECK_EQUAL(X[1], modulo(j + offset[1], N));
                BOOST_CHECK_EQUAL(X[2], modulo(k + offset[2], N));
            }
}
