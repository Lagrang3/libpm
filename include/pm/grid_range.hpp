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
        const data_fixed_t& _data;
        std::array<int, dim> _init, _end;
        int64_t num_elements = 1;
        size_t _size, _data_size;

       public:
        const_range(const data_fixed_t& data,
                    std::array<int, dim> init,
                    std::array<int, dim> end,
                    size_t size)
            : _data{data}, _init{init}, _end{end}, _size{size}
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
            : _data{that._data},
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

        int start(uint d) const
        {
            assert(d >= 0 and d < _init.size());
            return _init[d];
        }
        int stop(uint d) const
        {
            assert(d >= 0 and d < _end.size());
            return _end[d];
        }
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
        data_fixed_t& _data;
        std::array<int, dim> _init, _end;
        int64_t num_elements = 1;
        size_t _size, _data_size;

       public:
        range(data_fixed_t& data,
              std::array<int, dim> init,
              std::array<int, dim> end,
              size_t size)
            : _data{data}, _init{init}, _end{end}, _size{size}
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

        int start(uint d) const
        {
            assert(d >= 0 and d < _init.size());
            return _init[d];
        }
        int stop(uint d) const
        {
            assert(d >= 0 and d < _end.size());
            return _end[d];
        }
    };
}  // namespace PM
