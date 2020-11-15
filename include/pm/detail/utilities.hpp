#pragma once

#include <algorithm>
#include <cmath>

namespace PM
{
    namespace constants
    {
        const double pi = acos(-1.0);
    }

    namespace utilities
    {
        inline int decompose(const int Length, const int parts, const int rank)
        // Precondition: parts>0, Length>0, parts>rank>=0
        {
            int q = Length / parts, r = Length % parts;
            return q + (rank < r ? 1 : 0);
        }

        inline int sum_decompose(const int Length,
                                 const int parts,
                                 const int rank)
        // Precondition: parts>0, Length>0, parts>rank>=0
        {
            int q = Length / parts, r = Length % parts;
            return q * rank + std::min(r, rank);
        }

        inline int rank_decompose(const int Length,
                                  const int parts,
                                  const int index)
        // Precondition: parts>0, Length>0, Length>index>=0
        {
            int q = Length / parts, r = Length % parts;
            int Q = q + (r ? 1 : 0);
            int ans1 = index / Q, ans2 = (index - Q * r) / q + r;
            return ans1 < r ? ans1 : ans2;
        }

        template <int dim, class T>
        constexpr T power(T x)
        {
            T res = 1;
            for (int i = 0; i < dim; ++i)
                res *= x;
            return res;
        }

        inline int modulo(int x, int y)
        // Precondition: y>0
        {
            int r = x % y;
            return r < 0 ? r + std::abs(y) : r;
        };

        inline size_t index(int x,
                            int y,
                            int z,
                            const std::array<int64_t, 3>& strides) noexcept
        {
            return z * strides[2] + strides[1] * y + strides[0] * x;
        }

    }  // namespace utilities
}  // namespace PM
