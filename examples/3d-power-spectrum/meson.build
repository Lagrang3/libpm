three_d_power_spectrum_files= files(['main.cpp'])  
examples_files += three_d_power_spectrum_files

omp=dependency('openmp',required: false)
tbb=dependency('tbb', required: false)
pstl=dependency('pstl', required: false)
fftw3_omp = dependency('fftw3_omp', required: false)

if (omp.found() and tbb.found() and pstl.found() and fftw3_omp.found())
    executable('3d-pw',three_d_power_spectrum_files,
        include_directories: [libpm_include], 
        dependencies: [fftw3_omp,omp,tbb,pstl],
        link_with: libpm)
endif
