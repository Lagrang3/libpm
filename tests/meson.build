mpirun = find_program('mpirun')

test_files += files('utils.cpp')
utils = executable('utils','utils.cpp',
    include_directories: libpm_include,
    dependencies: boost_ut,
    link_with: libpm)
test('Utilities',utils)

test_files += files('fftw.cpp')
fftw = executable('fftw','fftw.cpp',
    include_directories: libpm_include,
    dependencies: [boost_ut,fftw3],
    link_with: libpm)
test('FFTW',fftw)
    
test_files += files('field.cpp')
field = executable('field',['field.cpp'],
    include_directories: libpm_include,
    dependencies: [boost,boost_ut,boost_mpi,mpi,fftw3])
test('field API',mpirun,
   args: ['-np','6',field.full_path()],
   is_parallel: false)

test_files += files('filters.cpp')
filters_ut = executable('filters',
    ['filters.cpp'],
    include_directories: libpm_include, 
    dependencies: [boost_ut],
    link_with: [libpm])
test('Filters',filters_ut)

# test_files += files('1d-interpolation.cpp')
# oned_interpolation_ut = executable('1d-interpolation-ut',
#     ['1d-interpolation.cpp'],
#     include_directories: [libpm_include], 
#     dependencies: [boost,fftw3,boost_ut],
#     link_with: [libpm])
# test('1D Interpolation',oned_interpolation_ut)

# 
