"""
Setup file for package `dsharp_opac`.
"""
import setuptools # noqa
from numpy.distutils.core import Extension
import pathlib

PACKAGENAME = 'dsharp_opac'

# the directory where this setup.py resides
HERE = pathlib.Path(__file__).parent


if __name__ == "__main__":
    from numpy.distutils.core import setup
    setup(name=PACKAGENAME,
          description='python routines to calculate mie opacities',
          version='1.1.1rc0',
          long_description=(HERE / "README.md").read_text(),
          long_description_content_type='text/markdown',
          url='https://github.com/birnstiel/dsharp_opac',
          author='Til Birnstiel & DSHARP collaboration',
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
