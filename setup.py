"""
Setup file for package `lp_opac`.
"""
from setuptools import setup
import os

PACKAGENAME = 'lp_opac'


setup(name=PACKAGENAME,
      description='python routines to calculate mie opacities',
      version='0.0.1',
      long_description=open(os.path.join(
          os.path.dirname(__file__), 'README.md')).read(),
      url='https://github.com/seanandrews/p484',
      author='LP',
      author_email='til.birnstiel@lmu.de',
      license='GPLv3',
      packages=[PACKAGENAME + '.py'],
      include_package_data=True,
      install_requires=['scipy', 'numpy', 'matplotlib'],
      zip_safe=False
      )
