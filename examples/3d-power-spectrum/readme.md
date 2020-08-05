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

Time in seconds:

|interpolation order|Poweri4 interlacing (s)|Poweri4 no interlacing (s)|libpm (s)|
|---|---|---|---|
|1|2.484|1.964|2.181|
|2|3.062|2.065|2.133|
|3|5.214|2.734|2.827|
|3|9.071|3.812|4.117|
