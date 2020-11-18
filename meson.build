project('libpm','cpp', version: '0.0.1',
    default_options: ['optimization=3','buildtype=release'])

if get_option('buildtype')=='release'
    add_global_arguments('-Ofast',language: 'cpp')
    add_global_arguments('-DNDEBUG',language: 'cpp')
endif

cpp=meson.get_compiler('cpp')
format = find_program('clang-format',required: false)
boost=dependency('boost')
boost_mpi=dependency('boost',modules :['mpi','serialization'])
boost_ut=dependency('boost',modules :['unit_test_framework'])
mpi = dependency('mpi',language: 'cpp')
mpirun = find_program('mpirun')

pkg=import('pkgconfig')
fftw3 = dependency('fftw3')

libpm_headers = []
subdir('include')

libpm = library('libpm',libpm_headers,
    include_directories: libpm_include,
    dependencies: [fftw3,mpi,boost_mpi,boost])
    
install_headers(libpm_headers)

pkg.generate(libpm,
    description: 'C++ library for Particle-Mesh cosmological codes')

examples_files=[]
subdir('examples')

test_files=[]
subdir('tests')

if (format.found())
    run_target('format', 
    command: [format,'-i',examples_files,libpm_headers,test_files])
endif
