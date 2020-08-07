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
        const data_vec_t& _data;
        std::array<int, dim> _state;
        int64_t _count;
        size_t _size;

       public:
        const_iterator(const const_range& p,
                       std::array<int, dim> init,
                       int64_t c)
            : parent_range{p},
              _data{p._data},
              _state{init},
              _count{c},
              _size{p._size}
        {
        }

        auto& operator++()
        {
            ++_count;

            // this first iteration out the loop
            // is done intentionally to gain some
            // performance
            ++_state[0];

            if (_state[0] < parent_range._end[0])
                return *this;

            _state[0] = parent_range._init[0];

            for (uint d = 1; d < dim; ++d)
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

        const auto& operator*()
        {
            const auto N = _size;

            size_t index = 0;
            for (int d = dim - 1; d >= 0; --d)
            {
                index = index * N + modulo(_state[d], N);
            }
            assert(index >= 0 and index < _data.size());
            return _data[index];
        }
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
        data_vec_t& _data;
        std::array<int, dim> _state;
        int64_t _count;
        size_t _size;

       public:
        iterator(const range& p, std::array<int, dim> init, int64_t c)
            : parent_range{p},
              _data{p._data},
              _state{init},
              _count{c},
              _size{p._size}
        {
        }

        auto& operator++()
        {
            ++_count;

            // this first iteration out the loop
            // is done intentionally to gain some
            // performance
            ++_state[0];

            if (_state[0] < parent_range._end[0])
                return *this;

            _state[0] = parent_range._init[0];

            for (uint d = 1; d < dim; ++d)
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

        auto& operator*()
        {
            const auto N = _size;

            size_t index = 0;
            for (int d = dim - 1; d >= 0; --d)
            {
                index = index * N + modulo(_state[d], N);
            }
            assert(index >= 0 and index < _data.size());
            return _data[index];
        }
        int count(uint d) const { return _state[d] - parent_range._init[d]; }
    };
}  // namespace PM
