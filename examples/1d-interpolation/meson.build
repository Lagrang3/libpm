examples_files += files('main.cpp')

exe = executable('1d-interpolation','main.cpp',
    include_directories: [libpm_include], 
    dependencies: [boost,mpi,boost_mpi],
    link_with: libpm)

plot_script = files('plot.py')

subdir('assets')
