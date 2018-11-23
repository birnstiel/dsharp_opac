"""
Setup file for package `dsharp_opac`.
"""
from numpy.distutils.core import Extension, setup
import os

PACKAGENAME = 'dsharp_opac'


setup(name=PACKAGENAME,
      description='python routines to calculate mie opacities',
      version='1.0.0',
      long_description=open(os.path.join(
          os.path.dirname(__file__), 'README.md')).read(),
      url='https://github.com/birnstiel/dsharp_opac',
      author='Til Birnstiel',
      author_email='til.birnstiel@lmu.de',
      license='GPLv3',
      packages=[PACKAGENAME],
      include_package_data=True,
      install_requires=['scipy', 'numpy', 'matplotlib'],
      zip_safe=False,
      ext_modules=[
          Extension(name='dsharp_opac.bhmie_fortran', sources=['dsharp_opac/bhmie_fortran.f90']),
          Extension(name='dsharp_opac.fit_module', sources=['dsharp_opac/fit_module.f90']),
          ],
      )
