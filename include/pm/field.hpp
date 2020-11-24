#pragma once

#include <array>
#include <boost/mpi.hpp>
#include <boost/mpi/cartesian_communicator.hpp>
//#include <valarray>

#include <pm/detail/fft.hpp>
#include <pm/detail/utilities.hpp>
#include <pm/filters.hpp>

namespace PM
{
    enum class direction_t
    {
        X,
        Y
    };

    template <class T,
              class sampler_t = filters::CIC,
              class interpolator_t = filters::CIC>
    class Field
    {
        using communicator = boost::mpi::cartesian_communicator;

        boost::mpi::cartesian_communicator com;
        int64_t N_glob;
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

        sampler_t W_sampler;
        interpolator_t W_interpolator;

        void reduce_ghosts();
        void clear_ghosts();
        void local_sample(const std::array<double, 3> p, double mass);
        std::array<double, 3> to_local(std::array<double, 3> p) const;

        int process(std::array<double, 3> p) const;

       public:
        Field(MPI_Comm raw_com, int64_t N, std::array<int, 2> proc);
        Field(boost::mpi::communicator boost_com,
              int64_t N,
              std::array<int, 2> proc);

        size_t size() const { return data.size(); }
        const auto& extents() const { return N_loc; }
        const auto& offset() const { return start_pos; }
        int64_t global_extents() const { return N_glob; }

        auto& operator()(int x, int y, int z);
        const auto& operator()(int x, int y, int z) const;

        // the `at` methods throw exceptions for out of range indexing
        auto& at(int x, int y, int z);
        const auto& at(int x, int y, int z) const;

        void fft(FFT_type type = FFT_type::forward);

        void tranpose_yz();
        void tranpose_xz();

        auto get_communicator() const { return com; }
        constexpr int ghost_thick() const
        {
            return std::max(sampler_t::width, interpolator_t::width) / 2;
        }

        T interpolate(double x, double y, double z) const;

        void clear();  // set to zero
        void update_ghosts();

        void sample(const std::vector<std::array<double, 3>> particles);

        T sum() const;

        std::vector<T> serialize() const;

        // TODO power spectrum
    };
}  // namespace PM

// TODO should i include here?
#include <pm/detail/field.hpp>
