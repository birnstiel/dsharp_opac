"""
Setup file for package `lp_opac`.
"""
from numpy.distutils.core import Extension, setup
import os

PACKAGENAME = 'lp_opac'


setup(name=PACKAGENAME,
      description='python routines to calculate mie opacities',
      version='0.0.2',
      long_description=open(os.path.join(
          os.path.dirname(__file__), 'README.md')).read(),
      url='https://github.com/seanandrews/p484',
      author='LP',
      author_email='til.birnstiel@lmu.de',
      license='GPLv3',
      packages=[PACKAGENAME],
      include_package_data=True,
      install_requires=['scipy', 'numpy', 'matplotlib'],
      zip_safe=False,
      ext_modules=[Extension(name='lp_opac.bhmie_draine', sources=[
              'lp_opac/bhmie_draine.f90'])],
      )
