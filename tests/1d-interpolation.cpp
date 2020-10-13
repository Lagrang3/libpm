#include <boost/test/unit_test.hpp>
#include <cmath>
#include <pm.hpp>

namespace ut = boost::unit_test;
const double pi = acos(-1.0);

template <int N, class filter_t, class callable>
void test(const int N_points, const double tolerance, callable F)
{
    // constexpr int N = k_max * 2 + 1;
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
            std::bind(&test<12, PM::NGP_filter, decltype(F)>, 10'000, 0.2, F),
            "NGP 12modes"));
        add(BOOST_TEST_CASE_NAME(
            std::bind(&test<100, PM::NGP_filter, decltype(F)>, 10'000, 0.02, F),
            "NGP 100modes"));
        add(BOOST_TEST_CASE_NAME(
            std::bind(&test<1000, PM::NGP_filter, decltype(F)>, 10'000, 0.002,
                      F),
            "NGP 1000modes"));

        add(BOOST_TEST_CASE_NAME(
            std::bind(&test<12, PM::CIC_filter, decltype(F)>, 10'000, 0.15, F),
            "CIC 12modes"));
        add(BOOST_TEST_CASE_NAME(
            std::bind(&test<100, PM::CIC_filter, decltype(F)>, 10'000, 0.015,
                      F),
            "CIC 100modes"));
        add(BOOST_TEST_CASE_NAME(
            std::bind(&test<1000, PM::CIC_filter, decltype(F)>, 10'000, 0.0015,
                      F),
            "CIC 1000modes"));

        add(BOOST_TEST_CASE_NAME(
            std::bind(&test<12, PM::TSC_filter, decltype(F)>, 10'000, 0.15, F),
            "TSC 12modes"));
        add(BOOST_TEST_CASE_NAME(
            std::bind(&test<100, PM::TSC_filter, decltype(F)>, 10'000, 0.015,
                      F),
            "TSC 100modes"));
        add(BOOST_TEST_CASE_NAME(
            std::bind(&test<1000, PM::TSC_filter, decltype(F)>, 10'000, 0.0015,
                      F),
            "TSC 1000modes"));

        add(BOOST_TEST_CASE_NAME(
            std::bind(&test<12, PM::PCS_filter, decltype(F)>, 10'000, 0.15, F),
            "PCS 12modes"));
        add(BOOST_TEST_CASE_NAME(
            std::bind(&test<100, PM::PCS_filter, decltype(F)>, 10'000, 0.015,
                      F),
            "PCS 100modes"));
        add(BOOST_TEST_CASE_NAME(
            std::bind(&test<1000, PM::PCS_filter, decltype(F)>, 10'000, 0.0015,
                      F),
            "PCS 1000modes"));

        add(BOOST_TEST_CASE_NAME(
            std::bind(&test<12, PM::Gaussian_filter, decltype(F)>, 10'000, 0.25,
                      F),
            "Gaussian 12modes"));
        add(BOOST_TEST_CASE_NAME(
            std::bind(&test<100, PM::Gaussian_filter, decltype(F)>, 10'000,
                      0.025, F),
            "Gaussian 100modes"));
        add(BOOST_TEST_CASE_NAME(
            std::bind(&test<1000, PM::Gaussian_filter, decltype(F)>, 10'000,
                      0.0025, F),
            "Gaussian 1000modes"));

        add(BOOST_TEST_CASE_NAME(
            std::bind(&test<13, PM::LowPass_filter<6, 13>, decltype(F)>, 10'000,
                      0.1, F),
            "LowPass 6/6modes"));
        add(BOOST_TEST_CASE_NAME(
            std::bind(&test<100, PM::LowPass_filter<6, 100>, decltype(F)>,
                      10'000, 0.1, F),
            "LowPass 6/100modes"));
        add(BOOST_TEST_CASE_NAME(
            std::bind(&test<1000, PM::LowPass_filter<6, 1000>, decltype(F)>,
                      10'000, 0.1, F),
            "LowPass 6/1000modes"));

        add(BOOST_TEST_CASE_NAME(
            std::bind(&test<13, PM::Grid_filter<13>, decltype(F)>, 10'000, 0.1,
                      F),
            "Sinc 6modes"));
        add(BOOST_TEST_CASE_NAME(
            std::bind(&test<100, PM::Grid_filter<100>, decltype(F)>, 10'000,
                      0.01, F),
            "Sinc 100modes"));
        add(BOOST_TEST_CASE_NAME(
            std::bind(&test<1000, PM::Grid_filter<1000>, decltype(F)>, 10'000,
                      0.001, F),
            "Sinc 1000modes"));
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
            std::bind(&test<100, PM::NGP_filter, decltype(F)>, 10'000, 0.35, F),
            "NGP 100modes"));
        add(BOOST_TEST_CASE_NAME(
            std::bind(&test<1000, PM::NGP_filter, decltype(F)>, 10'000, 0.035,
                      F),
            "NGP 1000modes"));

        add(BOOST_TEST_CASE_NAME(
            std::bind(&test<100, PM::CIC_filter, decltype(F)>, 10'000, 0.03, F),
            "CIC 100modes"));
        add(BOOST_TEST_CASE_NAME(
            std::bind(&test<1000, PM::CIC_filter, decltype(F)>, 10'000, 0.001,
                      F),
            "CIC 1000modes"));

        add(BOOST_TEST_CASE_NAME(
            std::bind(&test<100, PM::TSC_filter, decltype(F)>, 10'000, 0.03, F),
            "TSC 100modes"));
        add(BOOST_TEST_CASE_NAME(
            std::bind(&test<1000, PM::TSC_filter, decltype(F)>, 10'000, 0.001,
                      F),
            "TSC 1000modes"));

        add(BOOST_TEST_CASE_NAME(
            std::bind(&test<100, PM::PCS_filter, decltype(F)>, 10'000, 0.04, F),
            "PCS 100modes"));
        add(BOOST_TEST_CASE_NAME(
            std::bind(&test<1000, PM::PCS_filter, decltype(F)>, 10'000, 0.001,
                      F),
            "PCS 1000modes"));

        add(BOOST_TEST_CASE_NAME(
            std::bind(&test<100, PM::Gaussian_filter, decltype(F)>, 10'000, 0.2,
                      F),
            "Gaussian 100modes"));
        add(BOOST_TEST_CASE_NAME(
            std::bind(&test<1000, PM::Gaussian_filter, decltype(F)>, 10'000,
                      0.0015, F),
            "Gaussian 1000modes"));

        add(BOOST_TEST_CASE_NAME(
            std::bind(&test<13, PM::LowPass_filter<6, 13>, decltype(F)>, 10'000,
                      1e-8, F),
            "LowPass 6/6modes"));
        add(BOOST_TEST_CASE_NAME(
            std::bind(&test<100, PM::LowPass_filter<6, 100>, decltype(F)>,
                      10'000, 1e-8, F),
            "LowPass 6/100modes"));
        add(BOOST_TEST_CASE_NAME(
            std::bind(&test<1000, PM::LowPass_filter<6, 1000>, decltype(F)>,
                      10'000, 1e-8, F),
            "LowPass 6/1000modes"));

        add(BOOST_TEST_CASE_NAME(
            std::bind(&test<13, PM::Grid_filter<13>, decltype(F)>, 10'000, 1e-8,
                      F),
            "Sinc 6modes"));
        add(BOOST_TEST_CASE_NAME(
            std::bind(&test<100, PM::Grid_filter<100>, decltype(F)>, 10'000,
                      1e-8, F),
            "Sinc 100modes"));
        add(BOOST_TEST_CASE_NAME(
            std::bind(&test<1000, PM::Grid_filter<1000>, decltype(F)>, 10'000,
                      1e-8, F),
            "Sinc 1000modes"));
    }
};

ut::test_suite* init_unit_test_suite(int argc, char* argv[])
{
    ut::framework::master_test_suite().add(new triangle_test_suite);
    ut::framework::master_test_suite().add(new sine_test_suite);
    return 0;
}
