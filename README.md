# Mie-Opacity Library

This repository contains a module to carry out Mie opacity calculations based
on a range of astrophysically relevant optical datasets and to calculate the
default opacity tables used in the ALMA Large Program on Substructure in
Protoplanetary Disks (PI: S. Andrews).

This package was created by T. Birnstiel.

## Installation

got into the base directory (where this README.md is) and run

    pip install .

Note: some of the optical constants rely on data files that will be downloaded
upon first use of that particular class.


## Test

You can find a ipython notebook file under /tests that uses some of the
functionality of this package.
