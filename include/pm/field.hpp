#pragma once

#include <array>
#include <boost/mpi.hpp>
#include <boost/mpi/cartesian_communicator.hpp>
//#include <valarray>

#include <pm/detail/fft.hpp>
#include <pm/detail/utilities.hpp>

namespace PM
{
    enum class direction_t
    {
        X,
        Y
    };

    template <class T, int Nghost = 1>
    class Field
    {
        using communicator = boost::mpi::cartesian_communicator;

        boost::mpi::cartesian_communicator com;
        size_t N_glob;
        std::array<int64_t, 3> N_loc, start_pos, strides;
        std::array<int64_t, 3> gsizes_x, gsizes_y, gstrides_x, gstrides_y,
            gstart_x_up, gstart_y_up, gstart_x_down, gstart_y_down;

        std::vector<T> data;
        std::vector<T> ghost_x_up, ghost_x_down, ghost_y_up, ghost_y_down;

        inline size_t index(int x, int y, int z) const noexcept
        {
            return utilities::index(x, y, z, strides);
        }
        auto split_data(direction_t dir) const;
        void local_fft(FFT_type);

        void delegated_constructor();

       public:
        Field(MPI_Comm raw_com, size_t N, std::array<int, 2> proc);
        Field(boost::mpi::communicator boost_com,
              size_t N,
              std::array<int, 2> proc);

        size_t size() const { return data.size(); }
        const auto& extents() const { return N_loc; }
        const auto& offset() const { return start_pos; }

        auto& operator()(int x, int y, int z);
        const auto& operator()(int x, int y, int z) const;

        // the `at` methods throw exceptions for out of range indexing
        auto& at(int x, int y, int z);
        const auto& at(int x, int y, int z) const;

        void fft(FFT_type type = FFT_type::forward);

        void tranpose_yz();
        void tranpose_xz();

        auto get_communicator() const { return com; }
        constexpr int ghost_thick() const { return Nghost; }

        void update_ghosts();

        // TODO serialize
    };
}  // namespace PM

// TODO should i include here?
#include <pm/detail/field.hpp>
