
libgadget = dependency('gadget',required: false)
po = dependency('boost', modules: ['program_options'])

mini_body_files = files([
    'main.cpp',
    'integration.hpp',
    'parser.hpp',
    'state.hpp'])
    
examples_files += mini_body_files

if (libgadget.found() and po.found())
    minibody = executable('mini-body',mini_body_files,
        include_directories: [libpm_include], 
        dependencies: [libgadget,po,fftw3],
        link_with: libpm)

endif
