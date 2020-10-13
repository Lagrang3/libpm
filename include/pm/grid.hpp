#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <complex>
#include <fftw3.h>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <pm/utilities.hpp>

namespace PM
{
    template <int dim /* dimension = 1,2,3 */,
              class T /* precision */,
              class sampler_t /* sampling filter  */,
              class interpolator_t /* interpolating filter */>
    class grid
    {
        using data_t = std::complex<double>;
        using data_vec_t = std::vector<data_t>;
        using data_fixed_t = const std::unique_ptr<data_t[]>;

        class range;
        class iterator;

        class const_range;
        class const_iterator;

        const size_t _size, kN, _data_size;
        data_fixed_t _data;

        sampler_t W_in;
        interpolator_t W_out;

        fftw_plan plan;
        fftw_plan plan_inv;

        void compute_force(std::vector<T>& Force,
                           const std::vector<T>& Position,
                           const T dx) const
        {
            std::cerr << "Compute force ...\n";
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
            std::cerr << "Compute potential ...\n";
            const int64_t k_max = k_nyquist(), N = size();
            const double pi = acos(-1.0);
            auto index_to_k = [k_max, N](int64_t index) {
                return (index <= k_max ? index : index - N);
            };

            _data[0] = 0;  // const mode
            for (size_t i = 1; i < _data_size; ++i)
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

        auto get_index_range(const std::array<T, dim>& pos,
                             data_fixed_t& data,
                             const int int_width = 0,
                             const double width = 0)
        {
            const T dx = int_width > 0 ? width : size() * 0.5;
            const int dN = int_width > 0 ? int_width : size();
            const size_t N = size();
            std::array<int, dim> lower_corner, upper_corner;

            for (uint d = 0; d < dim; ++d)
            {
                auto xo = pos[d] * N;
                lower_corner[d] = std::ceil(xo - dx);
                upper_corner[d] = lower_corner[d] + dN;
            }
            return range(data, lower_corner, upper_corner, size());
        }
        auto get_index_range(const std::array<T, dim>& pos,
                             const data_fixed_t& data,
                             const int int_width = 0,
                             const double width = 0.) const
        {
            const T dx = int_width > 0 ? width : size() * 0.5;
            const int dN = int_width > 0 ? int_width : size();
            const size_t N = size();
            std::array<int, dim> lower_corner, upper_corner;

            for (uint d = 0; d < dim; ++d)
            {
                auto xo = pos[d] * N;
                lower_corner[d] = std::ceil(xo - dx);
                upper_corner[d] = lower_corner[d] + dN;
            }
            return const_range(data, lower_corner, upper_corner, size());
        }

        template <class filter_t>
        void get_weights(std::array<double, filter_t::int_width * dim>& Wval,
                         const_range index_range,
                         const std::array<T, dim>& pos) const
        {
            filter_t W;
            const size_t N = size();
            // std::array<std::vector<double>, dim> Wval;

            for (uint d = 0, c = 0; d < dim; ++d)
            {
                auto xo = pos[d] * N;
                for (int i = index_range.start(d); i < index_range.stop(d); ++i)
                {
                    assert(c >= 0 and c < Wval.size());
                    Wval[c++] = W(xo - i);
                }
            }
            // return Wval;
        }

        template <class filter_t>
        auto filter_fft() const
        {
            const int N = size();
            data_vec_t Wk(N);
            filter_t W;
            for (int i = -N / 2; i <= N / 2; ++i)
            {
                Wk[(i + N) % N] = W(i);
            }
#ifdef _OPENMP
            fftw_plan_with_nthreads(6);
#endif
            fftw_plan filter_plan =
                fftw_plan_dft_1d(N, reinterpret_cast<fftw_complex*>(&Wk[0]),
                                 reinterpret_cast<fftw_complex*>(&Wk[0]),
                                 FFTW_FORWARD, FFTW_ESTIMATE);
            fftw_execute(filter_plan);
            fftw_destroy_plan(filter_plan);
            return Wk;
        }

       public:
        grid() = delete;

        // constructor
        grid(size_t sz)
            : _size{sz},
              kN{(sz - 1) / 2},
              _data_size{PM::power<dim>(_size)},
              _data{new data_t[_data_size]}
        {
            auto N = size();

#ifdef _OPENMP
            fftw_plan_with_nthreads(6);
#endif
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
            return *range(_data, pos, pos, size()).begin();
        }
        const auto& at(std::array<int, dim> pos) const
        {
            return *const_range(_data, pos, pos, size()).begin();
        }

        size_t size() const { return _size; }
        int k_nyquist() const { return kN; }

        /* knowing the field in the grid, interpolate to any point in the box */
        double interpolate(const std::array<T, dim> pos) const noexcept
        {
            auto index_range = get_index_range(
                pos, _data, interpolator_t::int_width, interpolator_t::width);
            std::array<double, interpolator_t::int_width * dim> Wval;
            get_weights<interpolator_t>(Wval, index_range, pos);

            // std::cerr << "Interpolate at " << pos[0] << "\n";
            // std::cerr << "index range: " << index_range.start(0) << " - " <<
            // index_range.stop(0) << "\n";

            double answer = 0;
            for (auto i = index_range.begin(); i != index_range.end(); ++i)
            {
                double W = 1;
                for (uint d = 0; d < dim; ++d)
                {
                    uint idx = i.count(d) + d * interpolator_t::int_width;
                    assert(idx >= 0 and idx < Wval.size());
                    W *= Wval[idx];
                }
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
            std::cerr << "Force computation ...\n";
            if (Force.size() != Position.size())
                throw std::runtime_error(std::string{__PRETTY_FUNCTION__} +
                                         ": Force.size() != Position.size()");

            if (Force.size() % dim)
                throw std::runtime_error(
                    std::string{__PRETTY_FUNCTION__} +
                    ": vector<> Position is not a multiple of dim.");

            sample_density(Position);

            fft();

            compute_potential();

            fftw_execute(plan_inv);  // phi_k -> phi_x

            compute_force(Force, Position, 1.0 / size());
        }

        /*
            compute the fourier transform on the current grid.
        */
        void fft()
        {
            std::cerr << "Executing FFT ...\n";
            fftw_execute(plan);  // rho_x -> rho_k
        }

        /*
            samples a density field from
            a sum of dirac deltas
        */
        void sample_density(const std::vector<T>& Position)
        {
            std::cerr << "Sampling density ...\n";
            std::array<double, sampler_t::int_width * dim> Wval;
            for (size_t p = 0; p < Position.size(); p += dim)
            {
                std::array<T, dim> pos;
                std::copy(&Position[p], &Position[p + dim], pos.begin());

                auto index_range = get_index_range(
                    pos, _data, sampler_t::int_width, sampler_t::width);
                get_weights<sampler_t>(Wval, index_range, pos);

                for (auto i = index_range.begin(); i != index_range.end(); ++i)
                {
                    double W = 1;
                    for (uint d = 0; d < dim; ++d)
                    {
                        uint idx = d * sampler_t::int_width + i.count(d);
                        assert(idx >= 0 and idx < Wval.size());
                        W *= Wval[idx];
                    }
                    *i += data_t(W, 0);
                }
            }
        }

        /*
            compute the amplitude of the modes
            on the grid
        */
        auto get_modes() const
        {
            std::cerr << "Counting modes ...\n";

            const int k_max = size() / 2;
            std::vector<T> modes(k_max, 0);
            std::vector<int> count(k_max, 0);
            std::array<T, dim> center;

            std::fill(center.begin(), center.end(), 0);

            auto index_range = get_index_range(center, _data);

            for (auto i = index_range.begin(); i != index_range.end(); ++i)
            {
                double k = 0;
                for (int d = 0; d < dim; ++d)
                    k += i._state[d] * i._state[d];
                size_t idx = int(sqrt(k));
                if (idx < modes.size())
                {
                    modes[idx] += std::norm(*i);
                    ++count[idx];
                }
            }

            for (int i = 0; i < k_max; ++i)
                modes[i] /= (count[i] > 0 ? count[i] : 1);
            return modes;
        }

        /*
            Correction of Fourier modes
            due to sampling window function.
        */
        void sample_correction()
        {
            std::cerr << "Sample correction ...\n";
            std::array<T, dim> center;
            std::fill(center.begin(), center.end(), 0);
            auto index_range = get_index_range(center, _data);
            auto Wk = filter_fft<sampler_t>();
            const int N = size();
            for (auto i = index_range.begin(); i != index_range.end(); ++i)
            {
                data_t w = 1;
                for (int d = 0; d < dim; ++d)
                {
                    w *= Wk[modulo(i._state[d], N)];
                }
                *i /= w;
            }
        }

        ~grid()
        {
            fftw_destroy_plan(plan);
            fftw_destroy_plan(plan_inv);
        }
    };
}  // namespace PM

#include <pm/grid_iterator.hpp>
#include <pm/grid_range.hpp>
