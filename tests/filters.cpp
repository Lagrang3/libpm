#define BOOST_TEST_MODULE Filters
#include <boost/test/unit_test.hpp>
#include <cmath>
#include <pm/filters.hpp>
#include <random>
#include <vector>

namespace ut = boost::unit_test;

/* function with finite modes below-eq 5 */
auto good_signal = [](double x) {
    static const double pi = acos(-1.0);
    const int mod1 = 2, mod2 = 5;
    return sin(2 * pi * mod1 * x) + 2 * sin(2 * pi * mod2 * x);
};

/* function with infinte modes */
auto bad_signal = [](double x) { return x; };

auto random_vec(int N)
{
    std::vector<double> V(N);
    std::mt19937 rng(11);
    std::uniform_real_distribution<double> dis(-5, 5);
    for (auto& x : V)
        x = dis(rng);
    return V;
}

template <class F>
auto signal_vec(int N, F& fun)
{
    std::vector<double> V(N);
    for (int i = 0; i < N; ++i)
        V[i] = fun(double(i) / N);
    return V;
}

template <class filter_t>
struct interpolator
{
    interpolator(const std::vector<double>& V) : V{V} {}

    double operator()(double x) const
    //  0 <= x < 1
    {
        const int N = V.size();
        double r = 0;
        for (int j = 0; j < N; ++j)
            r += V[j] * W(x * N - j);
        return r;
    }

    filter_t W{};
    std::vector<double> V;
};

template <class F>
double difference_pts(const F& f, const std::vector<double>& V)
{
    double r = std::numeric_limits<double>::min();
    const int N = V.size();
    for (int i = 0; i < N; ++i)
    {
        r = std::max(r, std::abs(V[i] - f(double(i) / N)));
    }
    return r;
}

template <class F, class G>
double difference_functions(int N, const F& f, const G& g)
{
    double r = std::numeric_limits<double>::min();
    for (int i = 0; i < N; ++i)
    {
        r = std::max(r, std::abs(g(double(i) / N) - f(double(i) / N)));
    }
    return r;
}

BOOST_AUTO_TEST_CASE(Shannon_filters)
{
    double diff;
    const double eps = 1e-13;
    const double some = 0.01;
    {
        constexpr int N = 11;
        BOOST_TEST_CHECKPOINT("gridsize N odd, Nyquist function");
        const auto Vx = signal_vec(N, good_signal);
        interpolator<PM::filters::Sinc<N>> sinc(Vx);
        interpolator<PM::filters::LowPass<(N - 1) / 2, N>> low(Vx);

        diff = difference_pts(sinc, Vx);
        BOOST_CHECK_SMALL(diff, eps);

        diff = difference_pts(low, Vx);
        BOOST_CHECK_SMALL(diff, eps);

        diff = difference_functions(100'000, sinc, good_signal);
        BOOST_CHECK_SMALL(diff, eps);

        diff = difference_functions(100'000, low, good_signal);
        BOOST_CHECK_SMALL(diff, eps);
    }
    {
        constexpr int N = 11;
        BOOST_TEST_CHECKPOINT("gridsize N odd, non-Nyquist function");
        const auto Vx = signal_vec(N, bad_signal);
        interpolator<PM::filters::Sinc<N>> sinc(Vx);
        interpolator<PM::filters::LowPass<(N - 1) / 2, N>> low(Vx);

        diff = difference_pts(sinc, Vx);
        BOOST_CHECK_SMALL(diff, eps);

        diff = difference_pts(low, Vx);
        BOOST_CHECK_SMALL(diff, eps);

        diff = difference_functions(100'000, sinc, bad_signal);
        BOOST_CHECK_GE(diff, some);  // failure! because the input signal does
        // not satisfy Shannon-Nyquist theorem, it has infinite modes

        diff = difference_functions(100'000, low, bad_signal);
        BOOST_CHECK_GE(diff, some);  // failure! because the input signal does
        // not satisfy Shannon-Nyquist theorem, it has infinite modes
    }
    {
        constexpr int N = 12;
        BOOST_TEST_CHECKPOINT("gridsize N even, Nyquist function");
        const auto Vx = signal_vec(N, good_signal);
        interpolator<PM::filters::Sinc<N>> sinc(Vx);
        interpolator<PM::filters::LowPass<(N - 1) / 2, N>> low(Vx);

        diff = difference_pts(sinc, Vx);
        BOOST_CHECK_SMALL(diff, eps);

        diff = difference_pts(low, Vx);
        BOOST_CHECK_SMALL(diff, eps);

        diff = difference_functions(100'000, sinc, good_signal);
        BOOST_CHECK_GE(diff, some);  // failure! because Sinc is not an exact
        // interpolator of functions for N even

        diff = difference_functions(100'000, low, good_signal);
        BOOST_CHECK_SMALL(diff, eps);
    }
    {
        constexpr int N = 12;
        BOOST_TEST_CHECKPOINT("gridsize N even, non-Nyquist function");
        const auto Vx = signal_vec(N, bad_signal);
        interpolator<PM::filters::Sinc<N>> sinc(Vx);
        interpolator<PM::filters::LowPass<(N - 1) / 2, N>> low(Vx);

        diff = difference_pts(sinc, Vx);
        BOOST_CHECK_SMALL(diff, eps);

        diff = difference_pts(low, Vx);
        BOOST_CHECK_GE(diff,
                       some);  // failure! because the signal does not come
        // from a finite modes and LowPass does not interpolate exactly the grid
        // points for N even

        diff = difference_functions(100'000, sinc, bad_signal);
        BOOST_CHECK_GE(diff, some);  // failure! because the input signal does
        // not satisfy Shannon-Nyquist theorem, it has infinite modes

        diff = difference_functions(100'000, low, bad_signal);
        BOOST_CHECK_GE(diff, some);  // failure! because the input signal does
        // not satisfy Shannon-Nyquist theorem, it has infinite modes
    }
}
