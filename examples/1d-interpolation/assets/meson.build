data = custom_target('example 1d-interpolation data',
    build_by_default: true,
    output: [
        'sine_cic.dat',
        'sine_ngp.dat',
        'sine_tsc.dat',
        'sine_pcs.dat',
        'sine_gaussian.dat',
        'sine_exact.dat',
        'sine_points.dat',
        'triangle_cic.dat',
        'triangle_ngp.dat',
        'triangle_tsc.dat',
        'triangle_pcs.dat',
        'triangle_gaussian.dat',
        'triangle_exact.dat',
        'triangle_points.dat'
        ],
    command: [exe,'@OUTDIR@'])

custom_target('example 1d-interpolation plots',
    build_by_default: true,
    output: [
        'sine_cic.png',
        'sine_ngp.png',
        'sine_tsc.png',
        'sine_pcs.png',
        'sine_gaussian.png',
        'triangle_cic.png',
        'triangle_ngp.png',
        'triangle_tsc.png',
        'triangle_pcs.png',
        'triangle_gaussian.png'
        ],
    input: plot_script,
    depends: data,
    command: [python,'@INPUT@','@OUTDIR@'])

