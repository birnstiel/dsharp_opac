from IPython import get_ipython
ipython = get_ipython()
if ipython is not None:
    ipython.magic(u'matplotlib inline')

import matplotlib.pyplot as plt
import numpy as np
import os
import urllib
from astropy import constants as c

#from disklab import opacity
import dsharp_opac as opacity

import aux_functions as aux
aux.set_style()

if not os.path.isdir('figures'):
    os.mkdir('figures')
