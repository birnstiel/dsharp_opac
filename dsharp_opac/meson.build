python_sources = ['__init__.py', 'bhmie_python.py', 'dsharp_opac.py']
py3.install_sources(python_sources, subdir: 'dsharp_opac')

bhmie_fortran_sources = [
  'bhmie_fortran.f90'
  ]
bhmie_source = custom_target('bhmie_fortranmodule.c',
  input: bhmie_fortran_sources,
  output: ['bhmie_fortranmodule.c'],
  command: [py3, '-m', 'numpy.f2py', '--lower', '-m', 'bhmie_fortran', '@INPUT@', '--build-dir', 'dsharp_opac']
  )
py3.extension_module('bhmie_fortran',
  bhmie_fortran_sources+bhmie_source,
  incdir_f2py / 'fortranobject.c',
  include_directories: inc_np,
  dependencies: py_dep,
  install: true,
  subdir: 'dsharp_opac')

fit_module_fortran_sources = [
  'fit_module.f90'
  ]
fit_module_source = custom_target('fit_modulemodule.c',
  input: fit_module_fortran_sources,
  output: [
    'fit_modulemodule.c',
    'fit_module-f2pywrappers2.f90',
    ],
  command: [py3, '-m', 'numpy.f2py', '--lower', '-m', 'fit_module', '@INPUT@', '--build-dir', 'dsharp_opac']
  )
py3.extension_module('fit_module',
  fit_module_fortran_sources+fit_module_source,
  incdir_f2py / 'fortranobject.c',
  include_directories: inc_np,
  dependencies: py_dep,
  install: true,
  subdir: 'dsharp_opac')

subdir('data')
subdir('optical_constants')