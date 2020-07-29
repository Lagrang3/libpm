project('libmesh','cpp', version: '0.0.1',
    default_options: ['buildtype=debug'])

subdir('include')

_format = find_program('clang-format',required: false)
fftw3 = dependency('fftw3')
pkg=import('pkgconfig')

libpm = library('libpm',libpm_headers,include_directories: libpm_include)
install_headers(libpm_headers)

pkg.generate(libpm,description: 'C++ library for Particle-Mesh cosmological codes')

examples_files=[]
subdir('examples')

if (_format.found())
    run_target('format', command: [_format,'-i',examples_files,libpm_headers])
endif
