[![DOI](https://zenodo.org/badge/137751482.svg)](https://zenodo.org/badge/latestdoi/137751482)

# DSHARP Mie-Opacity Library

This repository contains a module to carry out Mie opacity calculations based on
a range of astrophysically relevant optical datasets and to calculate the
default opacity tables used in the Disk Substructures at High Angular Resolution
Project (DSHARP, PI: S. Andrews).

This package was created by T. Birnstiel and C.P. Dullemond.

## Credit

If you use this package, please cite Birnstiel et al. (2018).

If you use any optical constants, please cite the appropriate experimental papers. The reference will be displayed upon import of the constants and is also available as the

## Installation

-   From [source](https://github.com/birnstiel/dsharp_opac/archive/master.zip): go into the base directory (where this `README.md` is) and run

        pip install .

      Should you work on the sources, e.g. for including new optical constants or new mehtods, then we recommend linking the installation to the source folder. Changes to the sources will then be reflected upon `import dsharp_opac` without the need to reinstall the package.

-   Coming soon: from [pypi.org](pypi.org):

        pip install dsharp_opac

Note: some of the optical constants rely on data files that will be downloaded
upon first use of that particular class.

## Tests & Examples

You can find some jupyter notebooks in the [notebooks folder](notebooks/index.ipynb) that demonstrate some of the functionality of this package. It also contains the notebooks and data that were used to create the figures in Birnstiel et al. (2018).
