one_d_interpolation_files= files(['main.cpp'])  
examples_files += one_d_interpolation_files

one_d_interpolation = executable('1d-interpolation',one_d_interpolation_files,
    include_directories: [libpm_include], 
    dependencies: [fftw3],
    link_with: libpm)




