# DSHARP Mie-Opacity Library

This repository contains a module to carry out Mie opacity calculations based on
a range of astrophysically relevant optical datasets and to calculate the
default opacity tables used in the Disk Substructures at High Angular Resolution
Project (DSHARP, PI: S. Andrews).

This package was created by T. Birnstiel and C.P. Dullemond

## Installation

got into the base directory (where this README.md is) and run

    pip install .

Note: some of the optical constants rely on data files that will be downloaded
upon first use of that particular class.


## Test

You can find ipython notebooks under [/tests](/tests) that use some of the
functionality of this package as well as the IPython notebooks that were used to
create the figures in Birnstiel et al. (2018).

## Credit

If you use this package, please cite Birnstiel et al. (2018).

If you use any optical constants, please cite the appropriate experimental papers.
