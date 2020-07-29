
libgadget = dependency('gadget')

mini_body_files = files(['main.cpp','integration.hpp','parser.hpp','state.hpp'])  
examples_files += mini_body_files

minibody = executable('mini-body',mini_body_files,
    include_directories: [libpm_include], 
    dependencies: [libgadget,fftw3],
    link_with: libpm)



