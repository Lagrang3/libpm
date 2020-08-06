project('libpm','cpp', version: '0.0.1',
    default_options: ['optimization=3','buildtype=release'])

if get_option('buildtype')=='release'
    add_global_arguments('-Ofast',language: 'cpp')
    add_global_arguments('-DNDEBUG',language: 'cpp')
endif

subdir('include')

cpp=meson.get_compiler('cpp')
omp=dependency('openmp',required: false)
tbb=dependency('tbb', required: false)
pstl=dependency('pstl', required: false)

_format = find_program('clang-format',required: false)
fftw3 = dependency('fftw3_omp')
pkg=import('pkgconfig')

libpm = library('libpm',libpm_headers,
    include_directories: libpm_include,
    dependencies: [fftw3,omp,tbb,pstl])
    
install_headers(libpm_headers)

pkg.generate(libpm,description: 'C++ library for Particle-Mesh cosmological codes')

examples_files=[]
subdir('examples')

if (_format.found())
    run_target('format', command: [_format,'-i',examples_files,libpm_headers])
endif
