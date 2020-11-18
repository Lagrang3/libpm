#pragma once

#include <boost/mpi/collectives.hpp>
#include <boost/serialization/complex.hpp>
#include <boost/serialization/vector.hpp>
#include <cassert>
#include <pm/detail/fft.hpp>

namespace PM
{
    template <class T, class interpolator_t>
    void Field<T, interpolator_t>::delegated_constructor()
    {
        const auto coord = com.coordinates(com.rank());
        const auto top = com.topology();

        for (int i = 0; i < 2; ++i)
        {
            N_loc[i] = utilities::decompose(N_glob, top[i].size, coord[i]);
            start_pos[i] =
                utilities::sum_decompose(N_glob, top[i].size, coord[i]);
        }
        N_loc[2] = N_glob;
        start_pos[2] = 0;

        strides = {N_loc[1] * N_loc[2], N_loc[2], 1};

        data.resize(N_loc[0] * N_loc[1] * N_loc[2]);

        // ghost cells
        const int Nghost = ghost_thick();
        gsizes_x = {Nghost, 2 * Nghost + N_loc[1], N_loc[2]};
        gsizes_y = {N_loc[0], Nghost, N_loc[2]};

        gstrides_x = {gsizes_x[1] * gsizes_x[2], gsizes_x[2], 1};
        gstrides_y = {gsizes_y[1] * gsizes_y[2], gsizes_y[2], 1};

        ghost_x_down.resize(gsizes_x[0] * gsizes_x[1] * gsizes_x[2]);
        ghost_x_up.resize(ghost_x_down.size());

        ghost_y_down.resize(gsizes_y[0] * gsizes_y[1] * gsizes_y[2]);
        ghost_y_up.resize(ghost_y_down.size());

        gstart_x_down = {0 - Nghost, 0 - Nghost, 0};
        gstart_x_up = {N_loc[0], 0 - Nghost, 0};

        gstart_y_down = {0, 0 - Nghost, 0};
        gstart_y_up = {0, N_loc[1], 0};
    }

    template <class T, class interpolator_t>
    Field<T, interpolator_t>::Field(boost::mpi::communicator boost_com,
                                    size_t N,
                                    std::array<int, 2> proc)
        : com(boost_com,
              boost::mpi::cartesian_topology{{proc[0], true}, {proc[1], true}}),
          N_glob{N}
    {
        delegated_constructor();
    }

    template <class T, class interpolator_t>
    Field<T, interpolator_t>::Field(MPI_Comm raw_com,
                                    size_t N,
                                    std::array<int, 2> proc)
        : com(boost::mpi::communicator(
                  raw_com,
                  boost::mpi::comm_create_kind::comm_duplicate),
              boost::mpi::cartesian_topology{{proc[0], true}, {proc[1], true}}),
          N_glob{N}
    {
        delegated_constructor();
    }

    template <class T, class interpolator_t>
    void Field<T, interpolator_t>::fft(FFT_type type)
    {
        local_fft(type);

        tranpose_yz();
        local_fft(type);
        tranpose_yz();

        tranpose_xz();
        local_fft(type);
        tranpose_xz();
    }
    template <class T, class interpolator_t>
    void Field<T, interpolator_t>::local_fft(FFT_type type)
    {
        for (int i = 0; i < N_loc[0]; ++i)
            for (int j = 0; j < N_loc[1]; ++j)
                FFTW3(&data[index(i, j, 0)], &data[index(i, j + 1, 0)],
                      &data[index(i, j, 0)], type);
    }
    template <class T, class interpolator_t>
    auto Field<T, interpolator_t>::split_data(direction_t dir) const
    {
        const auto top = com.topology();
        const auto coord = com.coordinates(com.rank());
        const int com_size = top[dir == direction_t::X ? 0 : 1].size;

        std::vector<std::vector<T>> buff(com_size);
        for (int rank = 0; rank < com_size; ++rank)
        {
            const std::array<int64_t, 3> nloc{
                N_loc[0], N_loc[1],
                utilities::decompose(N_glob, com_size, rank)},
                strd{nloc[1] * nloc[2], nloc[2], 1};
            const auto start = utilities::sum_decompose(N_glob, com_size, rank);
            auto& V = buff[rank];

            V.resize(nloc[0] * nloc[1] * nloc[2]);

            for (int i = 0; i < nloc[0]; ++i)
                for (int j = 0; j < nloc[1]; ++j)
                    for (int k = 0; k < nloc[2]; ++k)
                        V[utilities::index(i, j, k, strd)] =
                            (*this)(i, j, k + start);
        }

        boost::mpi::communicator com_split;

        if (dir == direction_t::X)
            com_split = com.split(
                /*color=*/coord[1], /*rank=*/coord[0]);
        else
            com_split = com.split(
                /*color=*/coord[0], /*rank=*/coord[1]);

        boost::mpi::all_to_all(com_split, buff, buff);

        return buff;
    }

    template <class T, class interpolator_t>
    void Field<T, interpolator_t>::tranpose_yz()
    {
        const auto top = com.topology();
        const int com_size = top[1].size;
        const auto buff = split_data(direction_t::Y);

        for (int i = 0; i < N_loc[0]; ++i)
            for (int j = 0; j < N_loc[1]; ++j)
                for (int rank = 0; rank < com_size; ++rank)
                {
                    const std::array<int64_t, 3> nloc{
                        N_loc[0], utilities::decompose(N_glob, com_size, rank),
                        N_loc[1]};
                    const std::array<int64_t, 3> strd{nloc[2] * nloc[1],
                                                      nloc[2], 1};
                    const int start =
                        utilities::sum_decompose(N_glob, com_size, rank);
                    const auto& V = buff[rank];

                    for (int k = 0; k < nloc[1]; ++k)
                        (*this)(i, j, k + start) =
                            V[utilities::index(i, k, j, strd)];
                }
    }
    template <class T, class interpolator_t>
    void Field<T, interpolator_t>::tranpose_xz()
    {
        const auto top = com.topology();
        const int com_size = top[0].size;
        const auto buff = split_data(direction_t::X);

        for (int i = 0; i < N_loc[0]; ++i)
            for (int j = 0; j < N_loc[1]; ++j)
                for (int rank = 0; rank < com_size; ++rank)
                {
                    const std::array<int64_t, 3> nloc{
                        utilities::decompose(N_glob, com_size, rank), N_loc[1],
                        N_loc[0]};
                    const std::array<int64_t, 3> strd{nloc[2] * nloc[1],
                                                      nloc[2], 1};
                    const int start =
                        utilities::sum_decompose(N_glob, com_size, rank);
                    const auto& V = buff[rank];

                    for (int k = 0; k < nloc[0]; ++k)
                        (*this)(i, j, k + start) =
                            V[utilities::index(k, j, i, strd)];
                }
    }

    template <class T, class interpolator_t>
    auto& Field<T, interpolator_t>::operator()(int x, int y, int z)
    {
        // x planes comes first, because they are wider
        // +---------> y dir
        // | xxxxxx
        // | y    y
        // | y    y
        // | y    y
        // | xxxxxx
        // |
        // \/
        // x dir
        using utilities::index;
        z = utilities::modulo(z, N_glob);
        if (x < 0)
            return ghost_x_down[index(x, y, z, gstrides_x, gstart_x_down)];
        if (x >= N_loc[0])
            return ghost_x_up[index(x, y, z, gstrides_x, gstart_x_up)];
        if (y < 0)
            return ghost_y_down[index(x, y, z, gstrides_y, gstart_y_down)];
        if (y >= N_loc[1])
            return ghost_y_up[index(x, y, z, gstrides_y, gstart_y_up)];
        return data[index(x, y, z, strides)];
    }

    template <class T, class interpolator_t>
    const auto& Field<T, interpolator_t>::operator()(int x, int y, int z) const
    {
        using utilities::index;
        z = utilities::modulo(z, N_glob);
        if (x < 0)
            return ghost_x_down[index(x, y, z, gstrides_x, gstart_x_down)];
        if (x >= N_loc[0])
            return ghost_x_up[index(x, y, z, gstrides_x, gstart_x_up)];
        if (y < 0)
            return ghost_y_down[index(x, y, z, gstrides_y, gstart_y_down)];
        if (y >= N_loc[1])
            return ghost_y_up[index(x, y, z, gstrides_y, gstart_y_up)];
        return data[index(x, y, z, strides)];
    }
    template <class T, class interpolator_t>
    auto& Field<T, interpolator_t>::at(int x, int y, int z)
    {
        using utilities::index;
        z = utilities::modulo(z, N_glob);
        if (x < 0)
            return ghost_x_down.at(index(x, y, z, gstrides_x, gstart_x_down));
        if (x >= N_loc[0])
            return ghost_x_up.at(index(x, y, z, gstrides_x, gstart_x_up));
        if (y < 0)
            return ghost_y_down.at(index(x, y, z, gstrides_y, gstart_y_down));
        if (y >= N_loc[1])
            return ghost_y_up.at(index(x, y, z, gstrides_y, gstart_y_up));
        return data.at(index(x, y, z, strides));
    }
    template <class T, class interpolator_t>
    const auto& Field<T, interpolator_t>::at(int x, int y, int z) const
    {
        using utilities::index;
        z = utilities::modulo(z, N_glob);
        if (x < 0)
            return ghost_x_down.at(index(x, y, z, gstrides_x, gstart_x_down));
        if (x >= N_loc[0])
            return ghost_x_up.at(index(x, y, z, gstrides_x, gstart_x_up));
        if (y < 0)
            return ghost_y_down.at(index(x, y, z, gstrides_y, gstart_y_down));
        if (y >= N_loc[1])
            return ghost_y_up.at(index(x, y, z, gstrides_y, gstart_y_up));
        return data.at(index(x, y, z, strides));
    }

    template <class T, class interpolator_t>
    void Field<T, interpolator_t>::update_ghosts()
    {
        const int Nghost = ghost_thick();
        std::vector<T> buff;
        // exchange along Y, the X is wider
        auto next = com.shifted_ranks(1, 1);

        // high to low
        buff.resize(ghost_y_up.size());
        for (int i = 0, pos = 0; i < N_loc[0]; ++i)
            for (int j = 0; j < Nghost; ++j)
                for (int k = 0; k < N_loc[2]; ++k, ++pos)
                    buff[pos] = (*this)(i, j, k);

        com.sendrecv(next.first, 0x01, buff, next.second, 0x01, ghost_y_up);

        // low to high
        buff.resize(ghost_y_down.size());
        for (int i = 0, pos = 0; i < N_loc[0]; ++i)
            for (int j = N_loc[1] - Nghost; j < N_loc[1]; ++j)
                for (int k = 0; k < N_loc[2]; ++k, ++pos)
                    buff[pos] = (*this)(i, j, k);

        com.sendrecv(next.second, 0x00, buff, next.first, 0x00, ghost_y_down);

        // exchange along X
        next = com.shifted_ranks(0, 1);

        // high to low
        buff.resize(ghost_x_up.size());
        for (int i = 0, pos = 0; i < Nghost; ++i)
            for (int j = 0 - Nghost; j < N_loc[1] + Nghost; ++j)
                for (int k = 0; k < N_loc[2]; ++k, ++pos)
                    buff[pos] = (*this)(i, j, k);

        com.sendrecv(next.first, 0x11, buff, next.second, 0x11, ghost_x_up);

        // low to high
        buff.resize(ghost_x_down.size());
        for (int i = N_loc[0] - Nghost, pos = 0; i < N_loc[0]; ++i)
            for (int j = 0 - Nghost; j < N_loc[1] + Nghost; ++j)
                for (int k = 0; k < N_loc[2]; ++k, ++pos)
                    buff[pos] = (*this)(i, j, k);

        com.sendrecv(next.second, 0x10, buff, next.first, 0x10, ghost_x_down);
    }

    template <class T, class interpolator_t>
    T Field<T, interpolator_t>::interpolate(double x, double y, double z) const
    // x,y,z from 0 to N_glob
    {
        constexpr int width = interpolator_t::width;
        const double dx = 0.5 * width;
        std::array<double, 3 * width> Wval;
        const std::array<double, 3> coord{x, y, z};
        std::array<int, 3> low_coord;

        for (int d = 0; d < 3; ++d)
        {
            auto beg = Wval.begin() + width * d;
            low_coord[d] = std::ceil(coord[d] - dx);
            filters::weights(W_inter, low_coord[d], coord[d], beg, beg + width);
        }
        double ret = 0;

        for (int i = 0; i < width; ++i)
            for (int j = 0; j < width; ++j)
                for (int k = 0; k < width; ++k)
                {
                    double w = (*this)(i + low_coord[0], j + low_coord[1],
                                       k + low_coord[2]);
                    w *= Wval[i];
                    w *= Wval[j + width];
                    w *= Wval[k + 2 * width];
                    ret += w;
                }

        return ret;
    }
}  // namespace PM
