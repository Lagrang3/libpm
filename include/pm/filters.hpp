#pragma once
#include <cmath>
#include <pm/detail/utilities.hpp>

namespace PM
{
    namespace filters
    {
        class Gaussian
        {
            const double norm_factor = 1 / sqrt(2.0 * constants::pi);

           public:
            constexpr static int width = 8;

            double operator()(double x) const
            {
                return exp(-x * x * 0.5) * norm_factor;
            }
        };

        class NGP
        {
           public:
            constexpr static int width = 2;

            double operator()(double x) const
            {
                auto s = std::abs(x);
                return s < 0.5 ? 1 : (s > 0.5 ? 0 : 0.5);
            }
        };

        class CIC
        {
           public:
            constexpr static int width = 2;

            double operator()(double x) const
            {
                auto s = std::abs(x);
                return s < 1.0 ? 1 - s : 0;
            }
        };
        class TSC
        {
           public:
            constexpr static int width = 3;

            double operator()(double x) const
            {
                double s = std::abs(x);
                double r1 = 0.75 - s * s, r2 = (1.5 - s);
                return s < 0.5 ? r1 : (s < 1.5 ? r2 * r2 * 0.5 : 0);
            }
        };
        class PCS
        {
            const double norm_factor = 1. / 6;

           public:
            constexpr static int width = 4;

            double operator()(double x) const
            {
                double s = std::abs(x);
                double r1 = 4 - 6 * s * s + 3 * s * s * s, r2 = 2 - s;
                return norm_factor *
                       (s < 1.0 ? r1 : (s < 2.0 ? r2 * r2 * r2 : 0));
            }
        };

        template <int k_max, int N>
        class LowPass
        {
           public:
            constexpr static int width = N;

            double operator()(double x) const
            {
                const double iN = 1.0 / N;
                double s = iN * x * constants::pi;

                if (std::fabs(s) < 1e-8)
                    return iN * (2 * k_max + 1);

                double num = sin(s * (2 * k_max + 1)), den = sin(s);

                return num / den * iN;
            }
        };

        template <int N>
        class Sinc
        {
           public:
            constexpr static int width = N;

            double operator()(double x) const
            {
                const double iN = 1.0 / N;
                double s = iN * x * constants::pi;

                if (std::fabs(s) < 1e-8)
                    return 1.0;

                double num = sin(s * N), den = sin(s);

                return num / den * iN;
            }
        };

    }  // namespace filters
    // template <int N>
    // using Grid_filter = LowPass_filter<(N - 1) / 2, N>;

    // template <int N>
    // using Grid_filter = Sinc_filter< N>;

};  // namespace PM
