#include <boost/test/unit_test.hpp>
#include <cmath>
#include <pm.hpp>

namespace ut = boost::unit_test;
const double pi = acos(-1.0);

template <int k_max, class filter_t, class callable>
void test(const int N_points, const double tolerance, callable F)
{
    constexpr int N = k_max * 2 + 1;
    PM::grid<1, double, PM::Gaussian_filter, filter_t> mygrid(N);

    for (int i = 0; i < N; ++i)
        mygrid.at({i}) = F(double(i) / N);

    double max_diff = std::numeric_limits<double>::min();
    for (int i = 0; i < N_points; ++i)
    {
        double x = double(i) / N_points;
        auto diff = std::abs(mygrid.interpolate({x}) - F(x));
        max_diff = std::max(max_diff, diff);
    }
    BOOST_CHECK_SMALL(max_diff, tolerance);
}

template <class filter_t, class callable>
void test_case(callable F)
{
    test<12, filter_t>(10'000, 0.1, F);
    test<100, filter_t>(10'000, 0.01, F);
    test<1000, filter_t>(10'000, 0.001, F);
}

struct triangle_test_suite : public ut::test_suite
{
    triangle_test_suite() : ut::test_suite("Triangular shape")
    {
        auto F = [](double x) { return 3 * std::min(x, 1 - x); };

        add(BOOST_TEST_CASE_NAME(
            std::bind(&test_case<PM::NGP_filter, decltype(F)>, F), "NGP"));
        add(BOOST_TEST_CASE_NAME(
            std::bind(&test_case<PM::CIC_filter, decltype(F)>, F), "CIC"));
        add(BOOST_TEST_CASE_NAME(
            std::bind(&test_case<PM::TSC_filter, decltype(F)>, F), "TSC"));
        add(BOOST_TEST_CASE_NAME(
            std::bind(&test_case<PM::PCS_filter, decltype(F)>, F), "PCS"));
        add(BOOST_TEST_CASE_NAME(
            std::bind(&test_case<PM::Gaussian_filter, decltype(F)>, F),
            "Gaussian"));
    }
};
struct sine_test_suite : public ut::test_suite
{
    sine_test_suite() : ut::test_suite("Sine shape")
    {
        auto F = [](double x) {
            return sin(x * 2 * pi * 6) + sin(x * 2 * pi * 5);
        };

        add(BOOST_TEST_CASE_NAME(
            std::bind(&test_case<PM::NGP_filter, decltype(F)>, F), "NGP"));
        add(BOOST_TEST_CASE_NAME(
            std::bind(&test_case<PM::CIC_filter, decltype(F)>, F), "CIC"));
        add(BOOST_TEST_CASE_NAME(
            std::bind(&test_case<PM::TSC_filter, decltype(F)>, F), "TSC"));
        add(BOOST_TEST_CASE_NAME(
            std::bind(&test_case<PM::PCS_filter, decltype(F)>, F), "PCS"));
        add(BOOST_TEST_CASE_NAME(
            std::bind(&test_case<PM::Gaussian_filter, decltype(F)>, F),
            "Gaussian"));
    }
};

ut::test_suite* init_unit_test_suite(int argc, char* argv[])
{
    ut::framework::master_test_suite().add(new triangle_test_suite);
    ut::framework::master_test_suite().add(new sine_test_suite);
    return 0;
}
