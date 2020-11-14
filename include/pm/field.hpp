#pragma once

#include <array>
#include <boost/mpi.hpp>
#include <valarray>

#include <pm/detail/fft.hpp>
#include <pm/detail/utilities.hpp>

namespace PM
{
    enum class direction_t
    {
        X,
        Y
    };

    template <class T>
    class Field
    {
        using communicator = boost::mpi::cartesian_communicator;

        boost::mpi::cartesian_communicator com;
        size_t N_glob;
        std::array<int64_t, 3> N_loc, start_pos, strides;
        std::valarray<T> data;

        // TODO member ghost cell

        inline size_t index(int x, int y, int z) const noexcept
        {
            return utilities::index(x, y, z, strides);
        }
        auto split_data(direction_t dir) const;
        void local_fft(FFT_type);

       public:
        Field(boost::mpi::cartesian_communicator com, size_t N);

        size_t size() const { return data.size(); }
        const auto& extents() const { return N_loc; }
        const auto& offset() const { return start_pos; }
        auto& operator()(int x, int y, int z) { return data[index(x, y, z)]; }
        const auto& operator()(int x, int y, int z) const
        {
            return data[index(x, y, z)];
        }

        void fft(FFT_type type = FFT_type::forward);

        void tranpose_yz();
        void tranpose_xz();

        // TODO get cell ref
        // TODO serialize
        // TODO ghost cell update
    };
}  // namespace PM

// TODO should i include here?
#include <pm/detail/field.hpp>
