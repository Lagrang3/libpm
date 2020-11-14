#pragma once

#include <boost/mpi/collectives.hpp>
#include <boost/serialization/complex.hpp>
#include <boost/serialization/vector.hpp>
#include <cassert>
#include <pm/detail/fft.hpp>

namespace PM
{
    template <class T>
    Field<T>::Field(boost::mpi::cartesian_communicator com, size_t N)
        : com{com}, N_glob{N}
    {
        const auto coord = com.coordinates(com.rank());
        const auto top = com.topology();

        for (int i = 0; i < 2; ++i)
        {
            N_loc[i] = utilities::decompose(N, top[i].size, coord[i]);
            start_pos[i] = utilities::sum_decompose(N, top[i].size, coord[i]);
        }
        N_loc[2] = N;
        start_pos[2] = 0;

        strides = {N_loc[1] * N_loc[2], N_loc[2], 1};

        data.resize(N_loc[0] * N_loc[1] * N_loc[2]);
    }
    // TODO
    template <class T>
    void Field<T>::fft(FFT_type type)
    {
        local_fft(type);

        tranpose_yz();
        local_fft(type);
        tranpose_yz();

        tranpose_xz();
        local_fft(type);
        tranpose_xz();
    }
    template <class T>
    void Field<T>::local_fft(FFT_type type)
    {
        for (int i = 0; i < N_loc[0]; ++i)
            for (int j = 0; j < N_loc[1]; ++j)
                FFTW3(&data[index(i, j, 0)], &data[index(i, j + 1, 0)],
                      &data[index(i, j, 0)], type);
    }
    template <class T>
    auto Field<T>::split_data(direction_t dir) const
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

    template <class T>
    void Field<T>::tranpose_yz()
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
    template <class T>
    void Field<T>::tranpose_xz()
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
}  // namespace PM
