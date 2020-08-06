Example: Power spectrum on 3-dimensions
=======================================

This example tests the basic interface of the 
grid class and the filters defined in libpm.
Particularly the sampling routines, computation of the FFT
and computation of the amplitude of the Fourier modes.
The input consists of a snapshot from a cosmological simmulation 
*z=0* with *M=256^3* particles in 3-dimensions.
A grid of size *N=256* is then used to compute the power spectrum.
The power spectrum obtained with the libpm
using the filters:
- Nearest-grid-point (NGP),
- Cloud-in-cell (CIC),
- Triangular-shaped-cloud (TSC),
- Piecewise-cubic-split (PCS),
- Gaussian,
is then compared with the power spectrum obtained
using Sefusatti's code "Poweri4".

![](./assets/power-3d.png)

Time to solution comparison, between Sefusattis code Poweri4 and
libpm on a 6 core Intel(R) Core(TM) i7-9750H CPU @ 2.60GHz.

|interpolation order|Poweri4 interlacing (s)|Poweri4 no interlacing (s)|libpm (s)|libpm tbb (s)|
|---|---|---|---|---|
|1|2.505|1.857|2.058|1.223|
|2|3.053|1.959|1.866|1.253|
|3|5.457|2.590|2.424|1.564|
|3|9.971|3.593|3.587|2.084|

