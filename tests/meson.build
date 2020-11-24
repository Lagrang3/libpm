
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
   args: ['-np','6',field],
   is_parallel: false)

test_files += files('filters.cpp')
filters_ut = executable('filters',
    ['filters.cpp'],
    include_directories: libpm_include, 
    dependencies: [boost_ut],
    link_with: [libpm])
test('Filters',filters_ut)

test_files += files('ghost.cpp')
ghost_ut = executable('ghost',['ghost.cpp'],
    include_directories: libpm_include,
    dependencies: [boost,boost_ut,boost_mpi,mpi,fftw3])
test('field ghost cells',mpirun,
   args: ['-np','6',ghost_ut],
   is_parallel: false)

test_files += files('interpolation.cpp')
interpolation_ut = executable('interpolation-ut',
    ['interpolation.cpp'],
    include_directories: [libpm_include], 
    dependencies: [boost,boost_ut,mpi,boost_mpi],
    link_with: [libpm])
test('Interpolation',mpirun,
   args: ['-np','2',interpolation_ut],
   is_parallel: false)

test_files += files('sampling.cpp')
sampling_ut = executable('sampling-ut',
    ['sampling.cpp'],
    include_directories: [libpm_include], 
    dependencies: [boost,boost_ut,mpi,boost_mpi],
    link_with: [libpm])
test('Sampling',mpirun,
   args: ['-np','4',sampling_ut],
   is_parallel: false)
