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
shape, and *N=20*.

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

![](./assets/triangle_gaussian.png)

Sine+sine
---------

In this case *F= sin(2 pi x)+3 sin(4 pi x)* and *N=20*.


Nearest grid point (ngp) function compared to 
the exact interpolation filter (Whittaker-Shannon formula)

![](./assets/sine_ngp.png)


Cloud-in-cell (cic) function compared to 
the exact interpolation filter (Whittaker-Shannon formula)

![](./assets/sine_cic.png)

Triangular-shaped-cloud (tsc) function compared to 
the exact interpolation filter (Whittaker-Shannon formula)

![](./assets/sine_tsc.png)

Piecewise-cubic-spline (pcs) function compared to 
the exact interpolation filter (Whittaker-Shannon formula)

![](./assets/sine_pcs.png)

Gaussian function compared to 
the exact interpolation filter (Whittaker-Shannon formula)

![](./assets/sine_gaussian.png)
