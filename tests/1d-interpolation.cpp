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

struct triangle_test_suite : public ut::test_suite
{
    triangle_test_suite() : ut::test_suite("Triangular shape")
    {
        auto F = [](double x) { return 3 * std::min(x, 1 - x); };
        add(BOOST_TEST_CASE_NAME(
            std::bind(&test<12, PM::NGP_filter, decltype(F)>, 10'000, 0.1, F),
            "NGP 12modes"));
        add(BOOST_TEST_CASE_NAME(
            std::bind(&test<100, PM::NGP_filter, decltype(F)>, 10'000, 0.01, F),
            "NGP 100modes"));
        add(BOOST_TEST_CASE_NAME(
            std::bind(&test<1000, PM::NGP_filter, decltype(F)>, 10'000, 0.001,
                      F),
            "NGP 1000modes"));

        add(BOOST_TEST_CASE_NAME(
            std::bind(&test<12, PM::CIC_filter, decltype(F)>, 10'000, 0.1, F),
            "CIC 12modes"));
        add(BOOST_TEST_CASE_NAME(
            std::bind(&test<100, PM::CIC_filter, decltype(F)>, 10'000, 0.01, F),
            "CIC 100modes"));
        add(BOOST_TEST_CASE_NAME(
            std::bind(&test<1000, PM::CIC_filter, decltype(F)>, 10'000, 0.001,
                      F),
            "CIC 1000modes"));

        add(BOOST_TEST_CASE_NAME(
            std::bind(&test<12, PM::TSC_filter, decltype(F)>, 10'000, 0.1, F),
            "TSC 12modes"));
        add(BOOST_TEST_CASE_NAME(
            std::bind(&test<100, PM::TSC_filter, decltype(F)>, 10'000, 0.01, F),
            "TSC 100modes"));
        add(BOOST_TEST_CASE_NAME(
            std::bind(&test<1000, PM::TSC_filter, decltype(F)>, 10'000, 0.001,
                      F),
            "TSC 1000modes"));

        add(BOOST_TEST_CASE_NAME(
            std::bind(&test<12, PM::PCS_filter, decltype(F)>, 10'000, 0.1, F),
            "PCS 12modes"));
        add(BOOST_TEST_CASE_NAME(
            std::bind(&test<100, PM::PCS_filter, decltype(F)>, 10'000, 0.01, F),
            "PCS 100modes"));
        add(BOOST_TEST_CASE_NAME(
            std::bind(&test<1000, PM::PCS_filter, decltype(F)>, 10'000, 0.001,
                      F),
            "PCS 1000modes"));

        add(BOOST_TEST_CASE_NAME(
            std::bind(&test<12, PM::Gaussian_filter, decltype(F)>, 10'000, 0.15,
                      F),
            "Gaussian 12modes"));
        add(BOOST_TEST_CASE_NAME(
            std::bind(&test<100, PM::Gaussian_filter, decltype(F)>, 10'000,
                      0.015, F),
            "Gaussian 100modes"));
        add(BOOST_TEST_CASE_NAME(
            std::bind(&test<1000, PM::Gaussian_filter, decltype(F)>, 10'000,
                      0.0015, F),
            "Gaussian 1000modes"));
    }
};
struct sine_test_suite : public ut::test_suite
{
    sine_test_suite() : ut::test_suite("Sine shape")
    {
        auto F = [](double x) {
            return sin(2 * pi * x * 6) + sin(2 * pi * x * 5);
        };
        add(BOOST_TEST_CASE_NAME(
            std::bind(&test<100, PM::NGP_filter, decltype(F)>, 10'000, 0.2, F),
            "NGP 100modes"));
        add(BOOST_TEST_CASE_NAME(
            std::bind(&test<1000, PM::NGP_filter, decltype(F)>, 10'000, 0.02,
                      F),
            "NGP 1000modes"));

        add(BOOST_TEST_CASE_NAME(
            std::bind(&test<100, PM::CIC_filter, decltype(F)>, 10'000, 0.01, F),
            "CIC 100modes"));
        add(BOOST_TEST_CASE_NAME(
            std::bind(&test<1000, PM::CIC_filter, decltype(F)>, 10'000, 0.001,
                      F),
            "CIC 1000modes"));

        add(BOOST_TEST_CASE_NAME(
            std::bind(&test<100, PM::TSC_filter, decltype(F)>, 10'000, 0.01, F),
            "TSC 100modes"));
        add(BOOST_TEST_CASE_NAME(
            std::bind(&test<1000, PM::TSC_filter, decltype(F)>, 10'000, 0.001,
                      F),
            "TSC 1000modes"));

        add(BOOST_TEST_CASE_NAME(
            std::bind(&test<100, PM::PCS_filter, decltype(F)>, 10'000, 0.01, F),
            "PCS 100modes"));
        add(BOOST_TEST_CASE_NAME(
            std::bind(&test<1000, PM::PCS_filter, decltype(F)>, 10'000, 0.001,
                      F),
            "PCS 1000modes"));

        add(BOOST_TEST_CASE_NAME(
            std::bind(&test<100, PM::Gaussian_filter, decltype(F)>, 10'000,
                      0.03, F),
            "Gaussian 100modes"));
        add(BOOST_TEST_CASE_NAME(
            std::bind(&test<1000, PM::Gaussian_filter, decltype(F)>, 10'000,
                      0.0015, F),
            "Gaussian 1000modes"));
    }
};

ut::test_suite* init_unit_test_suite(int argc, char* argv[])
{
    ut::framework::master_test_suite().add(new triangle_test_suite);
    ut::framework::master_test_suite().add(new sine_test_suite);
    return 0;
}
