#define BOOST_TEST_MODULE Utilities
#include <boost/test/unit_test.hpp>
#include <pm/detail/utilities.hpp>

using namespace PM::utilities;

BOOST_AUTO_TEST_CASE(util_decompose)
{
    {
        const int L = 10, N = 4;
        BOOST_CHECK_EQUAL(rank_decompose(L, N, 0), 0);
        BOOST_CHECK_EQUAL(rank_decompose(L, N, 1), 0);
        BOOST_CHECK_EQUAL(rank_decompose(L, N, 2), 0);
        BOOST_CHECK_EQUAL(rank_decompose(L, N, 3), 1);
        BOOST_CHECK_EQUAL(rank_decompose(L, N, 4), 1);
        BOOST_CHECK_EQUAL(rank_decompose(L, N, 5), 1);
        BOOST_CHECK_EQUAL(rank_decompose(L, N, 6), 2);
        BOOST_CHECK_EQUAL(rank_decompose(L, N, 7), 2);
        BOOST_CHECK_EQUAL(rank_decompose(L, N, 8), 3);
        BOOST_CHECK_EQUAL(rank_decompose(L, N, 9), 3);
    }
    {
        const int L = 5, N = 7;
        BOOST_CHECK_EQUAL(decompose(L, N, 0), 1);
        BOOST_CHECK_EQUAL(decompose(L, N, 1), 1);
        BOOST_CHECK_EQUAL(decompose(L, N, 2), 1);
        BOOST_CHECK_EQUAL(decompose(L, N, 3), 1);
        BOOST_CHECK_EQUAL(decompose(L, N, 4), 1);
        BOOST_CHECK_EQUAL(decompose(L, N, 5), 0);
        BOOST_CHECK_EQUAL(decompose(L, N, 6), 0);

        BOOST_CHECK_EQUAL(sum_decompose(L, N, 0), 0);
        BOOST_CHECK_EQUAL(sum_decompose(L, N, 1), 1);
        BOOST_CHECK_EQUAL(sum_decompose(L, N, 2), 2);
        BOOST_CHECK_EQUAL(sum_decompose(L, N, 3), 3);
        BOOST_CHECK_EQUAL(sum_decompose(L, N, 4), 4);
        BOOST_CHECK_EQUAL(sum_decompose(L, N, 5), 5);
        BOOST_CHECK_EQUAL(sum_decompose(L, N, 6), 5);
    }
    {
        const int L = 5, N = 6;
        BOOST_CHECK_EQUAL(decompose(L, N, 0), 1);
        BOOST_CHECK_EQUAL(decompose(L, N, 1), 1);
        BOOST_CHECK_EQUAL(decompose(L, N, 2), 1);
        BOOST_CHECK_EQUAL(decompose(L, N, 3), 1);
        BOOST_CHECK_EQUAL(decompose(L, N, 4), 1);
        BOOST_CHECK_EQUAL(decompose(L, N, 5), 0);

        BOOST_CHECK_EQUAL(sum_decompose(L, N, 0), 0);
        BOOST_CHECK_EQUAL(sum_decompose(L, N, 1), 1);
        BOOST_CHECK_EQUAL(sum_decompose(L, N, 2), 2);
        BOOST_CHECK_EQUAL(sum_decompose(L, N, 3), 3);
        BOOST_CHECK_EQUAL(sum_decompose(L, N, 4), 4);
        BOOST_CHECK_EQUAL(sum_decompose(L, N, 5), 5);
    }
    {
        const int L = 5, N = 5;
        BOOST_CHECK_EQUAL(decompose(L, N, 0), 1);
        BOOST_CHECK_EQUAL(decompose(L, N, 1), 1);
        BOOST_CHECK_EQUAL(decompose(L, N, 2), 1);
        BOOST_CHECK_EQUAL(decompose(L, N, 3), 1);
        BOOST_CHECK_EQUAL(decompose(L, N, 4), 1);

        BOOST_CHECK_EQUAL(sum_decompose(L, N, 0), 0);
        BOOST_CHECK_EQUAL(sum_decompose(L, N, 1), 1);
        BOOST_CHECK_EQUAL(sum_decompose(L, N, 2), 2);
        BOOST_CHECK_EQUAL(sum_decompose(L, N, 3), 3);
        BOOST_CHECK_EQUAL(sum_decompose(L, N, 4), 4);
    }
    {
        const int L = 5, N = 4;
        BOOST_CHECK_EQUAL(decompose(L, N, 0), 2);
        BOOST_CHECK_EQUAL(decompose(L, N, 1), 1);
        BOOST_CHECK_EQUAL(decompose(L, N, 2), 1);
        BOOST_CHECK_EQUAL(decompose(L, N, 3), 1);

        BOOST_CHECK_EQUAL(sum_decompose(L, N, 0), 0);
        BOOST_CHECK_EQUAL(sum_decompose(L, N, 1), 2);
        BOOST_CHECK_EQUAL(sum_decompose(L, N, 2), 3);
        BOOST_CHECK_EQUAL(sum_decompose(L, N, 3), 4);
    }
    {
        const int L = 5, N = 3;
        BOOST_CHECK_EQUAL(decompose(L, N, 0), 2);
        BOOST_CHECK_EQUAL(decompose(L, N, 1), 2);
        BOOST_CHECK_EQUAL(decompose(L, N, 2), 1);

        BOOST_CHECK_EQUAL(sum_decompose(L, N, 0), 0);
        BOOST_CHECK_EQUAL(sum_decompose(L, N, 1), 2);
        BOOST_CHECK_EQUAL(sum_decompose(L, N, 2), 4);
    }
    {
        const int L = 5, N = 2;
        BOOST_CHECK_EQUAL(decompose(L, N, 0), 3);
        BOOST_CHECK_EQUAL(decompose(L, N, 1), 2);

        BOOST_CHECK_EQUAL(sum_decompose(L, N, 0), 0);
        BOOST_CHECK_EQUAL(sum_decompose(L, N, 1), 3);
    }
    {
        const int L = 5, N = 1;
        BOOST_CHECK_EQUAL(decompose(L, N, 0), 5);

        BOOST_CHECK_EQUAL(sum_decompose(L, N, 0), 0);
    }
}

BOOST_AUTO_TEST_CASE(util_modulo)
{
    {
        const int base = 1;
        BOOST_CHECK_EQUAL(modulo(-2, base), 0);
        BOOST_CHECK_EQUAL(modulo(-1, base), 0);
        BOOST_CHECK_EQUAL(modulo(0, base), 0);
        BOOST_CHECK_EQUAL(modulo(1, base), 0);
        BOOST_CHECK_EQUAL(modulo(2, base), 0);
    }

    {
        const int base = 2;
        BOOST_CHECK_EQUAL(modulo(-3, base), 1);
        BOOST_CHECK_EQUAL(modulo(-2, base), 0);
        BOOST_CHECK_EQUAL(modulo(-1, base), 1);
        BOOST_CHECK_EQUAL(modulo(0, base), 0);
        BOOST_CHECK_EQUAL(modulo(1, base), 1);
        BOOST_CHECK_EQUAL(modulo(2, base), 0);
        BOOST_CHECK_EQUAL(modulo(3, base), 1);
    }
    {
        const int base = 3;
        BOOST_CHECK_EQUAL(-4 % base, -1);
        BOOST_CHECK_EQUAL(modulo(-4, base), 2);
        BOOST_CHECK_EQUAL(modulo(-3, base), 0);
        BOOST_CHECK_EQUAL(modulo(-2, base), 1);
        BOOST_CHECK_EQUAL(modulo(-1, base), 2);
        BOOST_CHECK_EQUAL(modulo(0, base), 0);
        BOOST_CHECK_EQUAL(modulo(1, base), 1);
        BOOST_CHECK_EQUAL(modulo(2, base), 2);
        BOOST_CHECK_EQUAL(modulo(3, base), 0);
        BOOST_CHECK_EQUAL(modulo(4, base), 1);
    }
}

BOOST_AUTO_TEST_CASE(util_power)
{
    {
        const int base = 2;
        BOOST_CHECK_EQUAL(power<0>(base), 1);
        BOOST_CHECK_EQUAL(power<1>(base), 2);
        BOOST_CHECK_EQUAL(power<2>(base), 4);
        BOOST_CHECK_EQUAL(power<3>(base), 8);
        BOOST_CHECK_EQUAL(power<4>(base), 16);
        BOOST_CHECK_EQUAL(power<5>(base), 32);
    }
    {
        const int base = 3;
        BOOST_CHECK_EQUAL(power<0>(base), 1);
        BOOST_CHECK_EQUAL(power<1>(base), 3);
        BOOST_CHECK_EQUAL(power<2>(base), 9);
        BOOST_CHECK_EQUAL(power<3>(base), 27);
        BOOST_CHECK_EQUAL(power<4>(base), 81);
    }
    {
        const int base = -2;
        BOOST_CHECK_EQUAL(power<0>(base), 1);
        BOOST_CHECK_EQUAL(power<1>(base), -2);
        BOOST_CHECK_EQUAL(power<2>(base), 4);
        BOOST_CHECK_EQUAL(power<3>(base), -8);
        BOOST_CHECK_EQUAL(power<4>(base), 16);
    }
    {
        const int base = 1;
        BOOST_CHECK_EQUAL(power<0>(base), 1);
        BOOST_CHECK_EQUAL(power<1>(base), 1);
        BOOST_CHECK_EQUAL(power<2>(base), 1);
        BOOST_CHECK_EQUAL(power<3>(base), 1);
        BOOST_CHECK_EQUAL(power<4>(base), 1);
    }
    {
        const int base = -1;
        BOOST_CHECK_EQUAL(power<0>(base), 1);
        BOOST_CHECK_EQUAL(power<1>(base), -1);
        BOOST_CHECK_EQUAL(power<2>(base), 1);
        BOOST_CHECK_EQUAL(power<3>(base), -1);
        BOOST_CHECK_EQUAL(power<4>(base), 1);
    }
    {
        const int base = 0;
        BOOST_CHECK_EQUAL(power<0>(base), 1);
        BOOST_CHECK_EQUAL(power<1>(base), 0);
        BOOST_CHECK_EQUAL(power<2>(base), 0);
        BOOST_CHECK_EQUAL(power<3>(base), 0);
        BOOST_CHECK_EQUAL(power<4>(base), 0);
    }
}
