#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <complex>
#include <fftw3.h>
#include <stdexcept>
#include <string>
#include <vector>

namespace PM
{
    template <int dim /* dimension = 1,2,3 */,
              class T /* precision */,
              class sampler_t /* sampling filter  */,
              class interpolator_t /* interpolating filter */>
    class grid
    {
        class range;
        class iterator;

        class const_range;
        class const_iterator;

        size_t _size, kN;
        std::vector<std::complex<double>> _data;
        sampler_t W_in;
        interpolator_t W_out;

        fftw_plan plan;
        fftw_plan plan_inv;

        void compute_force(std::vector<T>& Force,
                           const std::vector<T>& Position,
                           const T dx) const
        {
            const T dx_inv = 1.0 / dx;
            for (size_t p = 0; p < Force.size(); p += dim)
            {
                std::array<T, dim> pos;
                std::copy(&Position[p], &Position[p + dim], pos.begin());

                for (size_t i = 0; i < dim; ++i)
                {
                    std::array<T, dim> pos_1{pos}, pos_2{pos};
                    pos_1[i] -= dx * 0.5;
                    pos_2[i] += dx * 0.5;
                    T phi_1 = interpolate(pos_1), phi_2 = interpolate(pos_2);

                    Force[p + i] = (phi_1 - phi_2) * dx_inv;
                }
            }
        }

        void compute_potential()
        {
            const int64_t k_max = k_nyquist(), N = size();
            auto index_to_k = [k_max, N](int64_t index) {
                return (index <= k_max ? index : index - N);
            };

            _data[0] = 0;  // const mode
            for (size_t i = 1; i < _data.size(); ++i)
            {
                int64_t k_squared = 0, tot = i;
                for (uint d = 0; d < dim; ++d)
                {
                    auto k = index_to_k(tot % N);
                    k_squared += k * k;
                    tot /= N;
                }
                _data[i] *= (-1.0) / (k_squared * pi);
            }
        }

        void sample_density(const std::vector<T>& Position)
        {
            for (size_t p = 0; p < Position.size(); p += dim)
            {
                std::array<T, dim> pos;
                std::copy(&Position[p], &Position[p + dim], pos.begin());

                auto index_range = get_index_range<sampler_t>(pos);
                auto Wval = get_weights<sampler_t>(index_range, pos);

                for (auto i = index_range.begin(); i != index_range.end(); ++i)
                {
                    double W = 1;
                    for (uint d = 0; d < dim; ++d)
                        W *= Wval[d][i.count(d)];

                    *i += std::complex<double>(W, 0);
                }
            }
        }

        template <class filter_t>
        auto get_index_range(const std::array<T, dim>& pos)
        {
            const T dx =
                filter_t::int_width > 0 ? filter_t::width : size() * 0.5;
            const int dN =
                filter_t::int_width > 0 ? filter_t::int_width : size();
            const size_t N = size();
            std::array<int, dim> lower_corner, upper_corner;

            for (uint d = 0; d < dim; ++d)
            {
                auto xo = pos[d] * N;
                lower_corner[d] = std::ceil(xo - dx);
                upper_corner[d] = lower_corner[d] + dN;
            }
            return grid<dim, T, sampler_t, interpolator_t>::range(
                *this, lower_corner, upper_corner);
        }
        template <class filter_t>
        auto get_index_range(const std::array<T, dim>& pos) const
        {
            const T dx =
                filter_t::int_width > 0 ? filter_t::width : size() * 0.5;
            const int dN =
                filter_t::int_width > 0 ? filter_t::int_width : size();
            const size_t N = size();
            std::array<int, dim> lower_corner, upper_corner;

            for (uint d = 0; d < dim; ++d)
            {
                auto xo = pos[d] * N;
                lower_corner[d] = std::ceil(xo - dx);
                upper_corner[d] = lower_corner[d] + dN;
            }
            return grid<dim, T, sampler_t, interpolator_t>::const_range(
                *this, lower_corner, upper_corner);
        }

        template <class filter_t>
        auto get_weights(const_range index_range,
                         const std::array<T, dim>& pos) const
        {
            filter_t W;
            const size_t N = size();
            std::array<std::vector<double>, dim> Wval;

            for (uint d = 0; d < dim; ++d)
            {
                auto xo = pos[d] * N;
                for (int i = index_range.start(d); i < index_range.stop(d); ++i)
                {
                    Wval[d].push_back(W(xo - i));
                }
            }
            return Wval;
        }

       public:
        grid() = delete;

        grid(size_t sz)
            : _size{sz}, kN{(sz - 1) / 2}, _data(_size * _size * _size)
        {
            auto N = size();

            switch (dim)
            {
                case 3:
                    plan = fftw_plan_dft_3d(
                        N, N, N, reinterpret_cast<fftw_complex*>(&_data[0]),
                        reinterpret_cast<fftw_complex*>(&_data[0]),
                        FFTW_FORWARD, FFTW_ESTIMATE);
                    plan_inv = fftw_plan_dft_3d(
                        N, N, N, reinterpret_cast<fftw_complex*>(&_data[0]),
                        reinterpret_cast<fftw_complex*>(&_data[0]),
                        FFTW_BACKWARD, FFTW_ESTIMATE);
                    break;
                case 2:
                    plan = fftw_plan_dft_2d(
                        N, N, reinterpret_cast<fftw_complex*>(&_data[0]),
                        reinterpret_cast<fftw_complex*>(&_data[0]),
                        FFTW_FORWARD, FFTW_ESTIMATE);
                    plan_inv = fftw_plan_dft_2d(
                        N, N, reinterpret_cast<fftw_complex*>(&_data[0]),
                        reinterpret_cast<fftw_complex*>(&_data[0]),
                        FFTW_BACKWARD, FFTW_ESTIMATE);
                    break;
                case 1:
                    plan = fftw_plan_dft_1d(
                        N, reinterpret_cast<fftw_complex*>(&_data[0]),
                        reinterpret_cast<fftw_complex*>(&_data[0]),
                        FFTW_FORWARD, FFTW_ESTIMATE);
                    plan_inv = fftw_plan_dft_1d(
                        N, reinterpret_cast<fftw_complex*>(&_data[0]),
                        reinterpret_cast<fftw_complex*>(&_data[0]),
                        FFTW_BACKWARD, FFTW_ESTIMATE);
                    break;
                default:

                    throw std::runtime_error(
                        "grid constructor Error: dim>3 is not supported");
            }
        }

        /* at function, accepts negative values in the input */
        auto& at(std::array<int, dim> pos)
        {
            auto modulo = [](int x, int y) {
                int r = x % y;
                return r < 0 ? r + std::abs(y) : r;
            };
            const auto N = size();

            size_t index = 0;
            for (int d = dim - 1; d >= 0; --d)
            {
                index = index * N + modulo(pos[d], N);
            }
            assert(index >= 0 and index < _data.size());
            return _data[index];
        }
        const auto& at(std::array<int, dim> pos) const
        {
            auto modulo = [](int x, int y) {
                int r = x % y;
                return r < 0 ? r + std::abs(y) : r;
            };
            const auto N = size();

            size_t index = 0;
            for (int d = dim - 1; d >= 0; --d)
            {
                index = index * N + modulo(pos[d], N);
            }
            assert(index >= 0 and index < _data.size());
            return _data[index];
        }

        size_t size() const { return _size; }
        int k_nyquist() const { return kN; }

        /* knowing the field in the grid, interpolate to any point in the box */
        double interpolate(const std::array<T, dim> pos) const
        {
            auto index_range = get_index_range<interpolator_t>(pos);
            auto Wval = get_weights<interpolator_t>(index_range, pos);

            // std::cerr << "Interpolate at " << pos[0] << "\n";
            // std::cerr << "index range: " << index_range.start(0) << " - " <<
            // index_range.stop(0) << "\n";

            double answer = 0;
            for (auto i = index_range.begin(); i != index_range.end(); ++i)
            {
                double W = 1;
                for (uint d = 0; d < dim; ++d)
                    W *= Wval[d][i.count(d)];
                // double dist = i.count(0) + index_range.start(0)- pos[0] *
                // size(); W = interpolator_t{}( dist/size() );
                answer += W * (*i).real();

                // std::cerr << " weight: " << W << '\n';
                // std::cerr << " distan: " << dist << '\n';
            }

            // std::cerr << " answer: " << answer << '\n';

            return answer;
        }

        /*
            this function is used to compute the Newtonian gravitational
            forces from the positions of the particles.
            (Gm = 1)
        */
        void operator()(std::vector<T>& Force, const std::vector<T>& Position)
        {
            if (Force.size() != Position.size())
                throw std::runtime_error(std::string{__PRETTY_FUNCTION__} +
                                         ": Force.size() != Position.size()");

            if (Force.size() % dim)
                throw std::runtime_error(
                    std::string{__PRETTY_FUNCTION__} +
                    ": vector<> Position is not a multiple of dim.");

            sample_density(Position);

            fftw_execute(plan);  // rho_x -> rho_k

            compute_potential();

            fftw_execute(plan_inv);  // phi_k -> phi_x

            compute_force(Force, Position, 1.0 / size());
        }

        ~grid()
        {
            fftw_destroy_plan(plan);
            fftw_destroy_plan(plan_inv);
        }
    };
}  // namespace PM

#include "libpm_grid_iterator.hpp"
#include "libpm_grid_range.hpp"
