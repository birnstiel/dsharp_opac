# Mie-Opacity scripts

This repository contains my Mie opacity scripts and the helper and testing functions. More details will follow ...

## Installation

got into `/tbhmie` and run `make download` then go back to the base directory and run e.g. `pip install -e .`.


## Test

to see if things are working try this from python/IPython:

    from tbhmie.tests import test_opac
    test_opac.test_opac()
