boost_ut=dependency('boost', modules : ['unit_test_framework'], 
    required: false)

examples_files += files(['1d-interpolation.cpp','filters.cpp'])  

if ( boost_ut.found() )
    oned_interpolation_ut = executable('1d-interpolation-ut',
        ['1d-interpolation.cpp'],
        include_directories: [libpm_include], 
        dependencies: [boost_ut,fftw3],
        link_with: [libpm])

    test('1D Interpolation',oned_interpolation_ut)
    
    filters_ut = executable('filters',
        ['filters.cpp'],
        include_directories: [libpm_include], 
        dependencies: [boost_ut,fftw3],
        link_with: [libpm])

    test('Filters',filters_ut)
endif
