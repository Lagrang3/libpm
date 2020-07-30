one_d_power_spectrum_files= files(['main.cpp'])  
examples_files += one_d_power_spectrum_files

executable('1d-pw',one_d_power_spectrum_files,
    include_directories: [libpm_include], 
    dependencies: [fftw3],
    link_with: libpm)

