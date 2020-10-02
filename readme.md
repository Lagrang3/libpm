libParticleMesh
===========================

This is a Particle-Mesh library for cosmological simulations,
written in C++.

Installation
------------

This code can be compiled and installed with the autoconfigure
tool `meson`. It is advised to create a build directory.

```
cd $PATH_TO_BUILD_DIR
meson $PATH_TO_SOURCES
meson configure --prefix $PATH_TO_INSTALLATION
ninja
```

Requirements
------------

- [`meson >= 0.54.2`](https://mesonbuild.com/)
- [`FFTW3`](http://fftw.org/)
- [`libgadget`](https://github.com/Lagrang3/libgadget)

Todo
----

- make an example with power spectrum using the
snapshots taken during gevolution-gadget simulations
and compare to seffusati
- run the mini-body example
- parallelize with threads
- parallelize with mpi
- implement unit tests
