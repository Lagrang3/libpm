#pragma once

#include <complex>
#include <fftw3.h>

namespace PM
{
    enum class FFT_type
    {
        forward,
        backward
    };

    int FFTW_sign(FFT_type T)
    {
        switch (T)
        {
            case FFT_type::forward:
                return FFTW_FORWARD;
            case FFT_type::backward:
                return FFTW_BACKWARD;
        }
        return -1;
    }

    template <class cRAiterator, class RAiterator>
    void FFTW3(cRAiterator in_beg,
               cRAiterator in_end,
               RAiterator out_beg,
               FFT_type T = FFT_type::forward)
    {
        const int n = std::distance(in_beg, in_end);
        fftw_plan plan;
        fftw_complex* data;

        data = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * n);
        plan = fftw_plan_dft_1d(n, data, data, FFTW_sign(T), FFTW_ESTIMATE);

        for (int i = 0; in_beg != in_end; ++i, ++in_beg)
        {
            data[i][0] = in_beg->real();
            data[i][1] = in_beg->imag();
        }

        fftw_execute(plan);

        for (int i = 0; i < n; ++i)
        {
            *out_beg = std::complex<double>{data[i][0], data[i][1]};
            ++out_beg;
        }

        fftw_free(data);
        fftw_destroy_plan(plan);
    }

}  // namespace PM
