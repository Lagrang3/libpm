three_d_power_spectrum_files= files(['main.cpp'])  
examples_files += three_d_power_spectrum_files

executable('3d-pw',three_d_power_spectrum_files,
    include_directories: [libpm_include], 
    dependencies: [fftw3],
    link_with: libpm)

