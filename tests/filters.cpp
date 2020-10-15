/*
    Unit Test:
    ----------

*/
#include <boost/test/unit_test.hpp>
#include <cmath>
#include <pm/filters.hpp>
#include <random>
#include <vector>

namespace ut = boost::unit_test;

template <int N, class filter_t>
void test_grid_pts(const std::vector<double>& V, const double tolerance)
{
    double max_diff = std::numeric_limits<double>::min();

    auto fun = [V](int i) {
        filter_t W;
        double r = 0;
        for (int j = 0; j < N; ++j)
            r += V[j] * W(i - j);
        return r;
    };

    for (int i = 0; i < N; ++i)
    {
        max_diff = std::max(max_diff, std::abs(V[i] - fun(i)));
    }
    BOOST_CHECK_SMALL(max_diff, tolerance);
}

/*
    Test the LowPass and Sinc filters
*/

auto random_vec(const int N)
{
    std::vector<double> V(N);
    std::mt19937 rng(11);
    std::uniform_real_distribution<double> dis(-5, 5);
    for (auto& x : V)
        x = dis(rng);
    return V;
}

template <int N>
auto sine_vec()
{
    std::vector<double> V(N);
    const double pi = acos(-1.0);
    constexpr int mod1 = 2, mod2 = 5;
    constexpr int kN = (N - 1) / 2;

    static_assert(mod1 <= kN);
    static_assert(mod2 <= kN);

    for (int i = 0; i < N; ++i)
        V[i] = sin(2 * pi * mod1) + 2 * sin(2 * pi * mod2);
    return V;
}

struct Shannon_test_suite : public ut::test_suite
{
    Shannon_test_suite() : ut::test_suite("Shannon Filters")
    {
        {
            constexpr int N = 11;
            const auto V{random_vec(N)};
            add(BOOST_TEST_CASE_NAME(
                std::bind(&test_grid_pts<N, PM::Sinc_filter<N>>, V, 1e-14),
                "Sinc " + std::to_string(N)));

            add(BOOST_TEST_CASE_NAME(
                std::bind(&test_grid_pts<N, PM::LowPass_filter<(N - 1) / 2, N>>,
                          V, 1e-14),
                "LowPass " + std::to_string(N)));
        }
        {
            constexpr int N = 12;
            const auto V{random_vec(N)};
            add(BOOST_TEST_CASE_NAME(
                std::bind(&test_grid_pts<N, PM::Sinc_filter<N>>, V, 1e-14),
                "Sinc " + std::to_string(N)));

            // Fails because for N even LowPass is not a delta kronecker
            // function and identiy only holds for specific functions F.
            // add(BOOST_TEST_CASE_NAME(
            //     std::bind(&test_grid_pts<N, PM::LowPass_filter<(N-1)/2,N> >,
            //     V,1e-14), "LowPass "+std::to_string(N)));
        }
        {
            constexpr int N = 12;
            const auto V{sine_vec<N>()};
            add(BOOST_TEST_CASE_NAME(
                std::bind(&test_grid_pts<N, PM::Sinc_filter<N>>, V, 1e-14),
                "Sinc (2) " + std::to_string(N)));

            // this test doesn't fail because the input function satisfies
            // F = F * W
            add(BOOST_TEST_CASE_NAME(
                std::bind(&test_grid_pts<N, PM::LowPass_filter<(N - 1) / 2, N>>,
                          V, 1e-14),
                "LowPass (2) " + std::to_string(N)));
        }
    }
};

ut::test_suite* init_unit_test_suite(int argc, char* argv[])
{
    ut::framework::master_test_suite().add(new Shannon_test_suite);
    return 0;
}
