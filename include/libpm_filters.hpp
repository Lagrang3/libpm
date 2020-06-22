#pragma once
#include <cmath>

namespace PM
{
    const double pi = acos(-1.0);

    class Gaussian_filter
    {
        const double norm_factor = 1 / sqrt(2.0 * pi);

       public:
        constexpr static int int_width = 8;
        constexpr static double width = int_width / 2;

        double operator()(double x) const
        {
            return exp(-x * x * 0.5) * norm_factor;
        }
    };

    class NGP_filter
    {
       public:
        constexpr static int int_width = 1;
        constexpr static double width = int_width / 2;

        double operator()(double x) const { return std::abs(x) < 0.5 ? 1 : 0; }
    };

    class CIC_filter
    {
       public:
        constexpr static int int_width = 2;
        constexpr static double width = int_width / 2;

        double operator()(double x) const
        {
            auto s = std::abs(x);
            return s < 1.0 ? 1 - s : 0;
        }
    };
    class TSC_filter
    {
       public:
        constexpr static int int_width = 3;
        constexpr static double width = int_width / 2;

        double operator()(double x) const
        {
            double s = std::abs(x);
            double r1 = 0.75 - s * s, r2 = (1.5 - s);
            return s < 0.5 ? r1 : (s < 1.5 ? r2 * r2 * 0.5 : 0);
        }
    };
    class PCS_filter
    {
        const double norm_factor = 1. / 6;

       public:
        constexpr static int int_width = 4;
        constexpr static double width = int_width / 2;

        double operator()(double x) const
        {
            double s = std::abs(x);
            double r1 = 4 - 6 * s * s + 3 * s * s * s, r2 = 2 - s;
            return norm_factor * (s < 1.0 ? r1 : (s < 2.0 ? r2 * r2 * r2 : 0));
        }
    };
};  // namespace PM
