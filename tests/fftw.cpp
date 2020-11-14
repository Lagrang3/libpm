#define BOOST_TEST_MODULE FFTW
#include <boost/test/unit_test.hpp>
#include <cmath>
#include <complex>
#include <pm/detail/fft.hpp>
#include <vector>

using namespace PM;
using cd = std::complex<double>;

#define DEBUG(x) std::cout << #x " = " << x.real() << ' ' << x.imag() << '\n'

BOOST_AUTO_TEST_CASE(FFT)
{
    const double tolerance = 1e-12;
    const double pi = acos(-1.0);
    {
        std::vector<cd> in{1};
        FFTW3(in.begin(), in.end(), in.begin());
        BOOST_CHECK_SMALL(std::abs(in[0] - cd{1}), tolerance);
    }
    {
        std::vector<cd> in{1, 0};
        FFTW3(in.begin(), in.end(), in.begin());
        BOOST_CHECK_SMALL(std::abs(in[0] - cd{1}), tolerance);
        BOOST_CHECK_SMALL(std::abs(in[1] - cd{1}), tolerance);
    }
    {
        std::vector<cd> in{1, 1};
        FFTW3(in.begin(), in.end(), in.begin(), FFT_type::backward);
        BOOST_CHECK_SMALL(std::abs(in[0] - cd{2}), tolerance);
        BOOST_CHECK_SMALL(std::abs(in[1] - cd{0}), tolerance);
    }
    {
        std::vector<cd> in{1, 2, 3};
        const cd e{cos(2 * pi / 3), -sin(2 * pi / 3)}, e2 = e * e;
        FFTW3(in.begin(), in.end(), in.begin());

        BOOST_CHECK_SMALL(std::abs(in[0] - (1.0 + 2.0 + 3.0)), tolerance);
        BOOST_CHECK_SMALL(std::abs(in[1] - (1.0 + 2.0 * e + 3.0 * e2)),
                          tolerance);
        BOOST_CHECK_SMALL(std::abs(in[2] - (1.0 + 2.0 * e2 + 3.0 * e2 * e2)),
                          tolerance);
    }
}
