#pragma once

#include <stdexcept>

namespace PM
{
    template <int dim /* dimension = 1,2,3 */,
              class T /* precision */,
              class sampler_t /* sampling filter  */,
              class interpolator_t /* interpolating filter */>
    class grid<dim, T, sampler_t, interpolator_t>::const_range
    {
        friend grid::const_iterator;

       protected:
        const grid& parent_grid;
        std::array<int, dim> _init, _end;
        int64_t num_elements = 1;

       public:
        const_range(const grid& p,
                    std::array<int, dim> init,
                    std::array<int, dim> end)
            : parent_grid{p}, _init{init}, _end{end}
        {
            for (uint d = 0; d < dim; ++d)
            {
                if (_end[d] < _init[d])
                    throw std::runtime_error(
                        "grid::range constructor Error: end<init");
                num_elements *= (_end[d] - _init[d]);
            }
        }
        const_range(range that)
            : parent_grid{that.parent_grid},
              _init{that._init},
              _end{that._end},
              num_elements{that.num_elements}
        {
        }

        auto begin() const { return grid::const_iterator(*this, _init, 0); }
        auto end() const
        {
            return grid::const_iterator(*this, _end, num_elements);
        }

        int start(uint d) const { return _init[d]; }
        int stop(uint d) const { return _end[d]; }
    };

    template <int dim /* dimension = 1,2,3 */,
              class T /* precision */,
              class sampler_t /* sampling filter  */,
              class interpolator_t /* interpolating filter */>
    class grid<dim, T, sampler_t, interpolator_t>::range
    {
        friend grid::iterator;
        friend grid::const_range;

       protected:
        grid& parent_grid;
        std::array<int, dim> _init, _end;
        int64_t num_elements = 1;

       public:
        range(grid& p, std::array<int, dim> init, std::array<int, dim> end)
            : parent_grid{p}, _init{init}, _end{end}
        {
            for (uint d = 0; d < dim; ++d)
            {
                if (_end[d] < _init[d])
                    throw std::runtime_error(
                        "grid::range constructor Error: end<init");
                num_elements *= (_end[d] - _init[d]);
            }
        }

        auto begin() { return grid::iterator{*this, _init, 0}; }
        auto end() { return grid::iterator{*this, _end, num_elements}; }

        int start(uint d) const { return _init[d]; }
        int stop(uint d) const { return _end[d]; }
    };
}  // namespace PM
