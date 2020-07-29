Example: Interpolation on 1-dimension
=====================================

This example tests the basic interface of the 
grid class and the filters defined in libpm.
Particularly the interpolation routines.
A function *F* in 1-dimension in the interval [0,1]
is sampled by a grid of size *N*, then 
*M* points are requested from the grid
interpolator, with *M>N*.

Triangular shape
----------------

In this case *F=3 min(x,1-x)* which is a triangular
shape, and *N=11*.

Nearest grid point (ngp) function compared to 
the exact interpolation filter (Whittaker-Shannon formula)

![](./assets/triangle_ngp.png)


Cloud-in-cell (cic) function compared to 
the exact interpolation filter (Whittaker-Shannon formula)

![](./assets/triangle_cic.png)

Triangular-shaped-cloud (tsc) function compared to 
the exact interpolation filter (Whittaker-Shannon formula)

![](./assets/triangle_tsc.png)

Piecewise-cubic-spline (pcs) function compared to 
the exact interpolation filter (Whittaker-Shannon formula)

![](./assets/triangle_pcs.png)

Gaussian function compared to 
the exact interpolation filter (Whittaker-Shannon formula)

![](./assets/triangle_gauss.png)

Low-pass filter function with 5 modes compared to 
the exact interpolation filter (Whittaker-Shannon formula)

![](./assets/triangle_low_pass.png)


Sine+sine
---------

In this case *F= sin(12 pi x)+sin(10 pi x)* and *N=100*.


Nearest grid point (ngp) function compared to 
the exact interpolation filter (Whittaker-Shannon formula)

![](./assets/sine_5_6_ngp.png)


Cloud-in-cell (cic) function compared to 
the exact interpolation filter (Whittaker-Shannon formula)

![](./assets/sine_5_6_cic.png)

Triangular-shaped-cloud (tsc) function compared to 
the exact interpolation filter (Whittaker-Shannon formula)

![](./assets/sine_5_6_tsc.png)

Piecewise-cubic-spline (pcs) function compared to 
the exact interpolation filter (Whittaker-Shannon formula)

![](./assets/sine_5_6_pcs.png)

Gaussian function compared to 
the exact interpolation filter (Whittaker-Shannon formula)

![](./assets/sine_5_6_gauss.png)

Low-pass filter function with 5 modes compared to 
the exact interpolation filter (Whittaker-Shannon formula).
In this case the input function has up to 6 modes
therefore the interpolating function (with only 5 modes)
can only *see* the *sin(10 pi x)* component.

![](./assets/sine_5_6_low_pass.png)
