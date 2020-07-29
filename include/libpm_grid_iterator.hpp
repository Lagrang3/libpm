#pragma once

namespace PM
{
    template <int dim /* dimension = 1,2,3 */,
              class T /* precision */,
              class sampler_t /* sampling filter  */,
              class interpolator_t /* interpolating filter */>
    class grid<dim, T, sampler_t, interpolator_t>::const_iterator
    {
        friend class grid<dim, T, sampler_t, interpolator_t>;

       protected:
        const const_range& parent_range;
        const grid& parent_grid;
        std::array<int, dim> _state;
        int64_t _count;

       public:
        const_iterator(const const_range& p,
                       std::array<int, dim> init,
                       int64_t c)
            : parent_range{p},
              parent_grid{p.parent_grid},
              _state{init},
              _count{c}
        {
        }

        auto& operator++()
        {
            ++_count;
            for (uint d = 0; d <= dim; ++d)
            {
                ++_state[d];
                if (_state[d] < parent_range._end[d])
                    break;

                _state[d] = parent_range._init[d];
            }

            return *this;
        }
        auto operator++(int)
        {
            iterator r = *this;
            ++(*this);
            return r;
        }
        bool operator==(const_iterator that) const
        {
            return _count == that._count;
        }
        bool operator!=(const_iterator that) const
        {
            return not(*this == that);
        }

        const auto& operator*() { return parent_grid.at(_state); }
        int count(uint d) const { return _state[d] - parent_range._init[d]; }
    };

    template <int dim /* dimension = 1,2,3 */,
              class T /* precision */,
              class sampler_t /* sampling filter  */,
              class interpolator_t /* interpolating filter */>
    class grid<dim, T, sampler_t, interpolator_t>::iterator
    {
        friend class grid<dim, T, sampler_t, interpolator_t>;

       protected:
        const range& parent_range;
        grid& parent_grid;
        std::array<int, dim> _state;
        int64_t _count;

       public:
        iterator(const range& p, std::array<int, dim> init, int64_t c)
            : parent_range{p},
              parent_grid{p.parent_grid},
              _state{init},
              _count{c}
        {
        }

        auto& operator++()
        {
            ++_count;
            for (uint d = 0; d <= dim; ++d)
            {
                ++_state[d];
                if (_state[d] < parent_range._end[d])
                    break;

                _state[d] = parent_range._init[d];
            }

            return *this;
        }
        auto operator++(int)
        {
            iterator r = *this;
            ++(*this);
            return r;
        }
        bool operator==(iterator that) const { return _count == that._count; }
        bool operator!=(iterator that) const { return not(*this == that); }

        auto& operator*() { return parent_grid.at(_state); }
        int count(uint d) const { return _state[d] - parent_range._init[d]; }
    };
}  // namespace PM
