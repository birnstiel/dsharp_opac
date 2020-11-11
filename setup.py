"""
Setup file for package `dsharp_opac`.
"""
import setuptools  # noqa
from numpy.distutils.core import Extension
import pathlib
import warnings

PACKAGENAME = 'dsharp_opac'

# the directory where this setup.py resides
HERE = pathlib.Path(__file__).parent


def read_version():
    with (HERE / PACKAGENAME / '__init__.py').open() as fid:
        for line in fid:
            if line.startswith('__version__'):
                delim = '"' if '"' in line else "'"
                return line.split(delim)[1]
        else:
            raise RuntimeError("Unable to find version string.")


if __name__ == "__main__":
    from numpy.distutils.core import setup

    extensions = [
        Extension(name='dsharp_opac.bhmie_fortran', sources=['dsharp_opac/bhmie_fortran.f90']),
        Extension(name='dsharp_opac.fit_module', sources=['dsharp_opac/fit_module.f90']),
    ]

    def setup_function(extensions):
        setup(
            name=PACKAGENAME,
            description='python routines to calculate mie opacities',
            version=read_version(),
            long_description=(HERE / "README.md").read_text(),
            long_description_content_type='text/markdown',
            url='https://github.com/birnstiel/dsharp_opac',
            author='Til Birnstiel & DSHARP collaboration',
            author_email='til.birnstiel@lmu.de',
            license='GPLv3',
            packages=[PACKAGENAME],
            package_dir={PACKAGENAME: 'dsharp_opac'},
            package_data={PACKAGENAME: [
                'optical_constants/*/*.*',
                'optical_constants/*/*/*.*',
                'data/*.*',
            ]},
            include_package_data=True,
            install_requires=['scipy', 'numpy', 'matplotlib', 'astropy', 'mpmath'],
            zip_safe=False,
            ext_modules=extensions
        )
    try:
        setup_function(extensions)
    except BaseException:
        warnings.warn('Mie calculations will use python routines -- this will be much slower than compiled code.')
        setup_function([])
