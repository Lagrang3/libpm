#pragma once

namespace PM
{
    template <int dim, class T>
    constexpr T power(T x)
    {
        T res = 1;
        for (int i = 0; i < dim; ++i)
            res *= x;
        return res;
    }
    int modulo(int x, int y)
    {
        int r = x % y;
        return r < 0 ? r + std::abs(y) : r;
    };
}  // namespace PM
