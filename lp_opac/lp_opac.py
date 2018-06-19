#!/usr/bin/env python
"""
This module contains opacity scripts and all the helper and testing routines.

- dielectric functions are objects, see bhmie.diel_*
  the functions `diel_*.nk` return the optical properties
"""
import numpy as np
import os
import sys
import warnings
import pkg_resources

try:
    from .bhmie_fortran import bhmie_fortran as bhmie_function
except ImportError:
    warnings.warn('could not import compiled mie code - mie calculation will be slow')
    from .bhmie_python import _bhmie_p

    def bhmie_function(x, nk, nangles):
        """
        wrapper for the python version to be callable just like the fortran version.

        Arguments:
        ----------

        x : float
            size parameter 2 pi a / lambda

        nk : complex
            complex ref. index = n + i * k, e.g. `complex(1.,0.)`

        nangles : int
            number of angles between 0 and 180 degree for which to return S1 & S2

        Output:
        -------
        S1, S2, Qext, Qabs, Qsca, Qback, gsca

        S1, S2 : arrays
            the matrix elements as function of angle

        Qext, Qabs, Qsca, Qback : float
            the extinction, absorption, scattering, backscattering coefficients

        gsca : float
            Henyey-Greenstein asymmetry factor
        """
        theta = np.linspace(0., 180., nangles + 1)
        return _bhmie_p(x, nk, theta)

try:
    import pymiecoated

    def bhmie_pymiecoated(x, nk, n_angles):
        """
        Wrapper to the Mie code of `pymiecoated`

        x : float
            size parameter x = 2 pi a / lambda

        nk : complex
            complex ref. index = n + i * k, e.g. `complex(1.,0.)`

        nangles : int
            number of angles between 0 and 180 degree for which to return S1 & S2

        Output:
        -------
        S1, S2, Qext, Qabs, Qsca, Qback, gsca

        S1, S2 : arrays
            the matrix elements as function of angle

        Qext, Qabs, Qsca, Qback : float
            the extinction, absorption, scattering, backscattering coefficients

        gsca : float
            Henyey-Greenstein asymmetry factor
        """
        #
        # the other code
        #
        opac = pymiecoated.Mie(x=x, m=nk)
        Qext = opac.qext()
        Qabs = opac.qabs()
        Qsca = opac.qsca()
        Qback = opac.qratio()
        gsca = opac.asy()

        S1 = np.zeros(n_angles + 1, dtype=complex)
        S2 = np.zeros(n_angles + 1, dtype=complex)

        angles = np.linspace(0., 180., n_angles + 1)

        for i, angle in enumerate(angles):
            S1[i], S2[i] = opac.S12(np.cos(angle / 180 * np.pi))

        return S1, S2, Qext, Qabs, Qsca, Qback, gsca
except ImportError:
    pass


def progress_bar(perc, text=''):
    """
    Displays a progress bar.

    This is a very simple progress bar which displays the given
    percentage of completion on the command line, overwriting its
    previous output.

    Arguments:
    perc    The percentage of completion (float), should be
             between 0 and 100. Only 100.0 finishes with the
             word "Done!".
    text    Possible text for describing the running process.

    Example:
    >>> import time
    >>> for i in linspace(0,100,1000):
    >>>     progress_bar(i,text='Waiting')
    >>>     time.sleep(0.005)
    """
    if text != '':
        text = text + ' ... '
    if perc == 100.0:
        sys.stdout.write('\r' + text + 'Done!\n')
        sys.stdout.flush()
    else:
        sys.stdout.write('\r' + text + '%d %%' % round(perc))
        sys.stdout.flush()


def download(packagedir):
    """
    Downloads the optical constants files. It works by getting the data from
    a json file `links.json` in `packagedir`.
    """
    import json
    import os
    from urllib.request import urlretrieve

    if not os.path.isdir(packagedir):
        packagedir = os.path.dirname(packagedir)

    with open(os.path.join(packagedir, 'links.json')) as f:
        data = json.load(f)
        for material, link in data.items():

            filename = link.split('/')[-1]

            print(f'material: {material}, downloading {filename}: ... ', end='')
            try:
                urlretrieve(link, filename=os.path.join(packagedir, filename))
                print('Done!')
            except Exception as ex:
                print('Failed!')
                raise ex


class diel_const(object):
    """
    Abstract class for dielectric constants objects
    """
    #
    # define the attributes of the class
    #
    datafile = ""
    _l = None
    _n = None
    _k = None
    _ll = None
    _ln = None
    _lk = None
    headerinfo = None
    material_str = None
    _lmin = None
    _lmax = None
    extrapol = False
    _has_negative_n = False

    def __init__(self):
        """
        Initialization of data and such
        """
        raise NameError('Abstract class not supposed to be initializable')

    def nk(self, l):
        """
        Return the optical properties interpolated at wavelength l

        Arguments:
        ----------

        l : float or array
        :    wavelength in cm

        Output:
        -------
        n : float
        :    real part of optical property

        k : float
        :    imaginary part of optical property
        """
        if np.array(l, ndmin=1).min() < self._lmin or np.array(l, ndmin=1).max() > self._lmax:
            raise NameError('{}: wavelength {:g} outside data-range [{:g},{:g}]'.format(type(self).__name__, l, self._lmin, self._lmax))

        log_interp = True
        if self._has_negative_n:
            # estimate if log interpolation is ok to do
            n = np.interp(l, self._l, self._n)
            if n < 0:
                log_interp = False

        if log_interp:
            result = 10.**np.interp(np.log10(l), self._ll, self._ln), 10.**np.interp(np.log10(l), self._ll, self._lk)
        else:
            result = np.interp(l, self._l, self._n), 10.**np.interp(np.log10(l), self._ll, self._lk)

        return result

    def extrapolate_constants(self, lmin, lmax):
        """
        Extend the data by extrapolation to longer wavelengths. Will start
        fitting for extrapolation at lmin and then extend the data up to lmax.
        """
        #
        # extrapolate
        #
        self.extrapol = True
        if self._has_negative_n:
            warnings.warn('Extrapolation for negative n values can cause issues')

        l_ext = np.logspace(np.log10(self._l[-1]), np.log10(lmax), 10)
        from scipy.optimize import curve_fit

        def f(x, a, b, c):
            return a + b * x + c * x**2

        i_min = abs(self._l - lmin).argmin()
        #
        # extrapolate n
        #
        res = curve_fit(f, np.log10(self._l[i_min:]), np.log10(self._n[i_min:]), [np.log10(self._n[-1]), 1, 0])
        n_ext = 10.**f(np.log10(l_ext), *res[0])
        n_new = np.append(self._n, n_ext)
        #
        # extrapolate k
        #
        res = curve_fit(f, np.log10(self._l[i_min:]), np.log10(self._k[i_min:]), [np.log10(self._k[-1]), 1, 0])
        k_ext = 10.**f(np.log10(l_ext), *res[0])
        k_new = np.append(self._k, k_ext)
        l_new = np.append(self._l, l_ext) # noqa

        # update attributes

        self._l = l_new
        self._n = n_new
        self._k = k_new
        self._ll = np.log10(self._l)
        self._ln = np.log10(self._n)
        self._lk = np.log10(self._k)
        self._lmin = self._l.min()
        self._lmax = self._l.max()
        self.material_str += ' - extrapolated'


class diel_from_lnk_file(diel_const):
    """
    returns the dielectric constants from a given lnk file.
    The lnk data needs to be in 3 columns:
    lambda [microns], n, k

    Arguments:
    ----------
    datafile : str
    :    path of lnk-file

    Keywords:
    ---------
    headerlines : int
    :    number of lines in header
    """

    def __init__(self, datafile, headerlines=0):
        """
        Overwrite the initialization of the parent class
        """
        #
        # open file, read header and data
        #
        self.datafile = datafile
        self.material_str = 'Optical constants from %s' % datafile
        f = open(self.datafile)
        self.headerinfo = [f.readline() for i in range(headerlines)]
        data = np.loadtxt(f)
        #
        # assign wavelength and optical constants
        #
        self._l = data[:, 0] * 1e-4
        self._n = data[:, 1]
        self._k = data[:, 2]
        self._ll = np.log10(self._l)
        self._ln = np.log10(self._n)
        self._lk = np.log10(self._k)
        self._lmin = self._l.min()
        self._lmax = self._l.max()

        if any(self._n <= 0):
            self._has_negative_n = True


class diel_henning(diel_const):
    """
    Returns the dielectric constants for the various optical constants
    from Thomas Hennings website.

    Arguments:
    ----------

    species : str
    :    name of the species to read. Possible choices are
         iron,olivine,organics,orthopyroxene,troilite,waterice


    Keywords:
    ---------

    type : str
    :    only for olivine and orthopyroxene decide between
         'low'    = low iron abundance
         'normal' = normal iron abundance
         'high'   = high iron abundance

    new : bool
    :    use the old version of the constants, for this case,
         only normal iron abundance is avaliable

    """

    def __init__(self, species, iron_abundance='normal', new=True):
        """
        Overwrite the initialization of the parent class
        """
        #
        # set the path and do some safety checks
        #
        self.material_str = species[0].upper() + species[1:] + ' (Henning)'

        if iron_abundance not in ['normal', 'low', 'high']:
            raise NameError('unknown iron abundance')
        if species not in ['iron', 'olivine', 'organics', 'orthopyroxene', 'troilite', 'waterice']:
            raise NameError('unknown abundance species')
        if (iron_abundance != 'normal') and ((not new) or (species not in ['olivine', 'orthopyroxene'])):
            raise NameError('low & high iron only available for new olivine and new orthopyroxene')
        #
        # set the file name
        #
        if species in ['iron', 'organics', 'troilite', 'waterice']:
            fname = species + new * 'k'
        elif species == 'olivine':
            if iron_abundance == 'normal':
                fname = new * 'olmg70k' + (not new) * 'olivine'
            elif iron_abundance == 'low':
                fname = 'olivine100k'
            elif iron_abundance == 'high':
                fname = 'olmg60k'
        elif species == 'orthopyroxene':
            if iron_abundance == 'normal':
                fname = new * 'pyrmg70k' + (not new) * 'orthopyr'
            elif iron_abundance == 'low':
                fname = 'pyrmg100k'
            elif iron_abundance == 'high':
                fname = 'pyrmg60k'

        self.datafile = pkg_resources.resource_filename(__name__, os.path.join(
            'optical_constants', 'henning', 'new' * new + 'old' * (not new), fname + '.lnk'))
        if not os.path.isfile(self.datafile):
            self.download()
        #
        # read data
        #
        print('Reading opacities from %s' % fname)
        data = np.loadtxt(self.datafile)
        #
        # assign wavelength and optical constants
        #
        self._l = data[:, 0] * 1e-4
        self._n = data[:, 1]
        self._k = data[:, 2]
        self._ll = np.log10(self._l)
        self._ln = np.log10(self._n)
        self._lk = np.log10(self._k)
        self._lmin = self._l.min()
        self._lmax = self._l.max()

    def download(self):
        for i in ['old', 'new']:
            path = pkg_resources.resource_filename(__name__, os.path.join('optical_constants', 'henning', i))
            download(path)


class diel_draine2003_astrosil(diel_const):
    """
    Returns the dielectric constants for astronomical silicates from
    Draine 2003. The data comes from callindex.out_sil.D03`
    which was downloaded from

        ftp://ftp.astro.princeton.edu/draine/dust/diel/callindex.out_silD03

    on 2014-01-21--17:16 EDT
    """

    def __init__(self):
        """
        Overwrite the initialization of the parent class
        """
        #
        # open file, read header and data
        #
        self.material_str = 'Astronomical Silicates (Draine 2003)'
        self.datafile = pkg_resources.resource_filename(__name__, os.path.join(
            'optical_constants', 'draine', 'callindex.out_silD03'))
        if not os.path.isfile(self.datafile):
            download(os.path.dirname(self.datafile))
        f = open(self.datafile)
        self.headerinfo = [f.readline() for i in range(5)]
        data = np.loadtxt(f)
        #
        # assign wavelength and optical constants
        #
        self._l = data[-1::-1, 0] * 1e-4
        self._n = data[-1::-1, 3] + 1.
        self._k = data[-1::-1, 4]
        self._ll = np.log10(self._l)
        self._ln = np.log10(self._n)
        self._lk = np.log10(self._k)
        self._lmin = self._l.min()
        self._lmax = self._l.max()


class diel_WD2001_astrosil(diel_const):
    """
    Returns the dielectric constants for astronomical silicates from
    Weingartner & Draine 2001. The data comes from file `eps_suvSil`
    which was downloaded from

        ftp://ftp.astro.princeton.edu/draine/dust/diel/eps_suvSil

    on 2018-06-18--11:20 EDT
    """

    def __init__(self):
        """
        Overwrite the initialization of the parent class
        """
        #
        # open file, read header and data
        #
        self.material_str = 'Astronomical Silicates (Weingartner & Draine 2001)'
        self.datafile = pkg_resources.resource_filename(__name__, os.path.join(
            'optical_constants', 'draine', 'eps_suvSil'))
        if not os.path.isfile(self.datafile):
            download(os.path.dirname(self.datafile))
        f = open(self.datafile)
        self.headerinfo = [f.readline() for i in range(9)]
        data = np.loadtxt(f)
        #
        # assign wavelength and optical constants
        #
        self._l = data[-1::-1, 0] * 1e-4
        self._n = data[-1::-1, 3] + 1.
        self._k = data[-1::-1, 4]
        self._ll = np.log10(self._l)
        self._ln = np.log10(self._n)
        self._lk = np.log10(self._k)
        self._lmin = self._l.min()
        self._lmax = self._l.max()


class diel_dl84_astrosil(diel_const):
    """
    Returns the dielectric constants for astronomical silicates from
    Draine & Lee 1984 (and Laor & Draine 1993). The data comes from eps_Sil
    which was downloaded from

        ftp://ftp.astro.princeton.edu/draine/dust/diel/eps_Sil

    on 2018-06-18--11:35 EDT
    """

    def __init__(self):
        """
        Overwrite the initialization of the parent class
        """
        #
        # open file, read header and data
        #
        self.material_str = 'Astronomical Silicates (Draine & Lee 1984)'
        self.datafile = pkg_resources.resource_filename(
            __name__, os.path.join('optical_constants', 'draine', 'eps_Sil'))
        f = open(self.datafile)
        if not os.path.isfile(self.datafile):
            download(os.path.dirname(self.datafile))
        self.headerinfo = [f.readline() for i in range(6)]
        data = np.loadtxt(f)
        #
        # assign wavelength and optical constants
        #
        self._l = data[::-1, 0] * 1e-4
        self._n = data[::-1, 1] + 1.
        self._k = data[::-1, 2]
        if any(self._n <= 0):
            self._has_negative_n = True
        self._ll = np.log10(self._l)
        self._ln = np.log10(self._n)
        self._lk = np.log10(self._k)
        self._lmin = self._l.min()
        self._lmax = self._l.max()


class diel_vacuum(diel_const):
    """
    Returns the dielectric constants for vacuum
    """
    extrapol = False

    def __init__(self):
        """
        Overwrite the initialization of the parent class
        """
        self.material_str = 'Vacuum'

    def nk(self, l):
        """
        Get the n and k values at the given wavelength.
        """
        l = np.array(l, ndmin=1) # noqa
        return np.array([np.ones(len(l)), np.zeros(len(l))])


class diel_zubko_carbon(diel_const):
    """
    Returns the dielectric constants for carbon grains from Zubko et
    al. 1996 (the BE, ACH2 values). The data was OCRed by Til Birnstiel, no
    guarantee for correctness.

    Arguments:
    ----------

    sample : str
        which of the samples to use, implemented are: 'ACH2', 'BE'

    Keywords:
    ---------

    extrapol : bool
        whether or not to extrapolate beyond ~0.19 mm

    lmax : float
        upper limit for extrapolation range
    """

    def __init__(self, sample='ACH2', extrapol=False, lmin=4.765e-2, lmax=1.0):
        """
        Overwrite the initialization of the parent class
        """
        #
        # open files, read data
        #
        directory = os.path.join('optical_constants', 'zubko+1996')

        self.material_str = 'Carbonaceous Grains (Zubko et al. 1996)'
        self.datafile = directory

        E = np.loadtxt(pkg_resources.resource_filename(
            __name__, os.path.join(self.datafile, f'zubko_E_{sample}.txt')))[-1::-1]
        n = np.loadtxt(pkg_resources.resource_filename(
            __name__, os.path.join(self.datafile, f'zubko_n_{sample}.txt')))[-1::-1]
        k = np.loadtxt(pkg_resources.resource_filename(
            __name__, os.path.join(self.datafile, f'zubko_k_{sample}.txt')))[-1::-1]
        l = 0.00012398419292004205 / E  # E = h*c/lambda in CGS # noqa
        #
        # assign wavelength and optical constants
        #
        self._l = l
        self._n = n
        self._k = k
        self._ll = np.log10(self._l)
        self._ln = np.log10(self._n)
        self._lk = np.log10(self._k)
        self._lmin = self._l.min()
        self._lmax = self._l.max()

        if extrapol:
            self.extrapolate_constants(lmin, lmax)


class diel_warren(diel_const):
    """
    Returns the dielectric constants for water ice according
    to the data by Warren 1984 or the newer version

    the new, updated data was downloaded from

        http://www.atmos.washington.edu/ice_optical_constants/

    on Jan 23, 2014 and is based on Warren & Brandt (2008).

    The old data was taken from the journal website, where
    the coldest temperature column was used.

    Keywords:
    ---------

    new : bool
        use the old version of the constants if true, old if false

    """

    def __init__(self, new=True):
        """
        Overwrite the initialization of the parent class
        """
        #
        # set the path and do some safety checks
        #
        if new:
            self.material_str = 'Water Ice (Warren & Brandt 2008)'
        else:
            self.material_str = 'Water Ice (Warren 1984)'
        #
        # set the file name
        #
        if new:
            fname = 'IOP_2008_ASCIItable.dat'
        else:
            fname = 'warren_1984.txt'
        self.datafile = pkg_resources.resource_filename(
            __name__, os.path.join('optical_constants', 'warren', fname))
        if not os.path.isfile(self.datafile):
            download(self.datafile)
        #
        # read data and assign wavelength and optical constants
        #
        data = np.loadtxt(self.datafile)
        self._l = data[:, 0] * 1e-4
        self._n = data[:, 1]
        self._k = data[:, 2]
        #
        # assign derived attributes
        #
        self._ll = np.log10(self._l)
        self._ln = np.log10(self._n)
        self._lk = np.log10(self._k)
        self._lmin = self._l.min()
        self._lmax = self._l.max()


class diel_luca(diel_const):
    """
    Returns the dielectric constants from Luca Riccis files
    (cf. Ricci et al. 2010, A&A vol. 512, p. 15)

    Arguments:
    ---------
    species : string
    :    which species to return; can be silicate, ice, or carbon

    Keywords:
    ---------
    extrapol : bool
    :    whether or not to extrapolate beyond the given data up lmax (default: 10cm)

    lmax : float
    :    upper limit of extrapolation range (upper extrapolation is quadratic)

    lmin : float
    :    lower limit of extrapolation range (lower extrapolation is constant)

    Notes:
    ------
    Default is no extrapolation. If extrapol==True, the upper range will be extrapolated
    to 10 cm or the value given as lmax. Lower extrapolation (constant value) will only
    be done if a lower limit is explicitly given as lmin.

    """

    def __init__(self, species, extrapol=False, lmin=None, lmax=10.0):
        """
        Overwrite the initialization of the parent class
        """
        #
        # set the path and do some safety checks
        #
        self.material_str = ('Lucas ' + species).replace('ice', 'water ice')
        self.extrapol = extrapol
        #
        # set the file name
        #
        if species == 'silicate':
            fname = 'silicate.dat'
        elif species == 'ice':
            fname = 'H2Oice_Warren.dat'
        elif species == 'carbon':
            fname = 'aC_ACH2_Zubko.dat'
        else:
            raise NameError('%s got unknown species: %s' % (type(self).__name__, species))
        self.datafile = pkg_resources.resource_filename(
            __name__, os.path.join('optical_constants', 'luca', fname))
        #
        # read data and assign wavelength and optical constants
        #
        f = open(self.datafile)
        self.headerinfo = [f.readline() for i in range(2)]
        data = np.loadtxt(f)
        l = data[:, 0] * 1e-4 # noqa
        n = data[:, 1] + 1.0
        k = data[:, 2]
        #
        # extrapolate
        #
        if extrapol:
            #
            # UPPER EXTRAPOLATION: log-quadratic
            #
            le = np.logspace(np.log10(l[-1]), np.log10(lmax), 10)
            from scipy.optimize import curve_fit

            def f(x, a, b, c):
                return a + b * x + c * x**2

            i_min = abs(l - 1e-1).argmin()
            #
            # extrapolate n
            #
            res = curve_fit(f, np.log10(l[i_min:]), np.log10(n[i_min:]), [np.log10(n[-1]), 1, 0])
            ne = 10.**f(np.log10(le), *res[0])
            n = np.append(n, ne)
            #
            # extrapolate k
            #
            res = curve_fit(f, np.log10(l[i_min:]), np.log10(k[i_min:]), [np.log10(k[-1]), 1, 0])
            ke = 10.**f(np.log10(le), *res[0])
            k = np.append(k, ke)
            l = np.append(l, le) # noqa
            #
            # LOWER EXTRAPOLATION: constant
            #
            if lmin is not None:
                n = np.append(n[0], n)
                k = np.append(k[0], k)
                l = np.append(lmin, l) # noqa
        #
        # assign the attributes
        #
        self._l = l
        self._n = n
        self._k = k
        self._ll = np.log10(self._l)
        self._ln = np.log10(self._n)
        self._lk = np.log10(self._k)
        self._lmin = self._l.min()
        self._lmax = self._l.max()


class diel_mixed(diel_const):
    """
    This is a dielectric_constant class that mixes the various
    dielectric constants given their abundances.

    Arguments:
    ----------
    constants : list of objects of class diel_const
    :    all the materials that should be mixed

    abundances : array
    :    the volume fractions for each material

    rule : str
    :    the mixing rule. Possible choices are
         'Bruggeman'
    """

    def __init__(self, constants, abundances, rule='Bruggeman', extrapol=False):
        """
        Initialize mixed dielectric constants.

        Arguments
        ---------
        constants : list
            a list of optical constants (diel_const)

        abundances : list
            volume fractions of the different constants

        rule : str
            the mixing rule; 'Bruggeman' or 'Maxwell-Garnett'

        extrapol : bool
            whether or not to extrapolate constants beyond there defined interval
        """
        if rule not in ['Bruggeman', 'Maxwell-Garnett']:
            raise NameError('Unknown mixing rule: %s' % rule)

        self.material_str = '%s-Mix of %i species' % (rule, len(constants))
        self.constants = constants
        self.abundances = abundances
        self.rule = rule
        self.extrapol = extrapol

    def nk(self, l):
        """
        Returns the mixed optical constants at the given wavelength

        Arguments:
        ----------

        l : float or array
        :    wavelength in cm

        Output:
        -------
        n : float
        :    real part of mixed optical property

        k : float
        :    imaginary part of mixed optical property
        """
        from mpmath import findroot
        l_arr = np.array(l, ndmin=1)
        eps_mean = np.empty(np.shape(l_arr)).astype('complex')

        for i, l in enumerate(l_arr):
            #
            # calculate eps = (n - I*k)**2 for each material
            #
            eps = np.array([complex(*c.nk(l)).conjugate()**2 for c in self.constants])

            if self.rule.lower() == 'bruggeman':
                #
                # define the mixing rule and solve for the mixed value
                #
                def fct(x):
                    return sum(self.abundances * ((eps - x) / (eps + 2 * x)))
                eps_mean[i] = complex(findroot(fct, complex(0.5, 0.5)))
            elif self.rule.lower() == 'maxwell-garnett':
                #
                # e.g. kataoka et al. 2014, eq. 3
                #
                fj_gammaj = np.array(self.abundances) * 3. / (eps + 2.)
                eps_mean[i] = (fj_gammaj * eps).sum() / (fj_gammaj.sum())
            #
            # return n and k
            #
            eps_mean = np.sqrt(eps_mean)
            return np.array([eps_mean.real.squeeze(), -eps_mean.imag.squeeze()])


def powerlaw_N_of_a(a, a_max, q, rho_s):
    """
    Gives the size distribution N(a) ~ a**-q, normalized to 1 g total
    mass with upper cut-off a_max.

    Arguments:
    ----------

    a : array
        the size grid in cm

    a_max : float
        the upper size cut-off in cm

    q : float
        the negative exponent, like 3.5 for MRN

    rho_s : float
        the material density in g cm**-3

    Output:
    -------
    dist : array
        the size distribution N(a)

    """
    m = 4. * np.pi / 3. * rho_s * a**3
    dist = a**-q * (a <= a_max)
    dist = dist / np.trapz(dist * m, x=a)
    return dist


def gaussian_N_of_a(a, a_mean, sigma_a, rho_s):
    """
    Gives a gaussian size distribution
    N(a) ~ exp(-(a - a_mean)**2 / (2 * sigma**2)),
    normalized to 1 g total mass.

    Arguments:
    ----------

    a : array
        the size grid in cm

    a_max : float
        the upper size cut-off in cm

    q : float
        the negative exponent, like 3.5 for MRN

    rho_s : float
        the material density in g cm**-3

    Output:
    -------
    dist : array
        the size distribution N(a)

    """
    m = 4. * np.pi / 3. * rho_s * a**3
    dist = np.exp(-(a - a_mean)**2 / (2 * sigma_a**2))
    dist = dist / np.trapz(dist * m, x=a)
    return dist


def get_kappa_from_q(a, m, q_abs, q_sca):
    """
    Converts absorption and scattering coefficients [unitless] to absorption and
    scattering opacities [cm^2/g].

    Arguments:
    ----------

    q_abs, q_sca : arrays
        absorption and scattering coefficients, shape (n_sizes, n_wavelength)

    a, m: arrays
        particle size and particle mass arrays

    Output:
    -------

    kappa_abs, kappa_sca : arrays
        opacities in units of [cm^2/g]
    """
    n_lam = q_abs.shape[1]
    kappa_abs = q_abs * np.tile(np.pi * a**2 / m, [n_lam, 1]).transpose()
    kappa_sca = q_sca * np.tile(np.pi * a**2 / m, [n_lam, 1]).transpose()
    return kappa_abs, kappa_sca


def get_size_averaged_opacity(a, lam, n, rho_s, diel_const=None, q_abs=None,
                              q_sca=None, k_abs=None, k_sca=None):
    """
    Averages the opacity over the given size distribution.

    Arguments:
    ----------

    a : array
        the size grid in cm

    lam : array
        the wavelength grid in cm

    n : array
        the size distribution N(a)

    rho_s : float
        the material density of the dust grains


    Keywords:
    ---------

    One of the following three options has to be given:

    diel_const : object of class diel_const
        the dielectric constants to be used

    q_abs,q_sca : array
        if the opacity coefficients for all sizes and wavelength has already
        been calculated, then you can pass it along, otherwise
        it will be calculated on the fly.

    q_abs,q_sca : array
        if the opacity for all sizes and wavelength has already
        been calculated, then you can pass it along, otherwise
        it will be calculated on the fly.

    Output:
    -------
    kappa_abs,kappa_sca : array
        the opacity at each wavelength averaged over the size
        distribution and normalized per 1 g of dust.
    """
    assert (diel_const is not None), "diel_const needs to be "

    if (q_abs is None or q_sca is None) and (k_abs is None or k_sca is None) \
            and (diel_const is None):
        raise AssertionError('Either (diel_const) or (q_abs, q_sca) or (k_abs, k_sca) needed as input')
    #
    # some basic conversions
    #
    m = 4. * np.pi / 3. * rho_s * a**3
    sig = n * m * a
    sig = sig / sum(sig)
    #
    # calculate the opacities ...
    #
    if (q_abs is None or q_sca is None) and (k_abs is None or k_sca is None):
        q_abs, q_sca = get_mie_coefficients(a, lam, diel_const)
    elif (k_abs is None or k_sca is None):
        kappa_abs, kappa_sca = get_kappa_from_q(a, m, q_abs, q_sca)
    #
    # ... average them over the size distribution ...
    #
    kappa_abs_m = np.zeros(len(lam))
    kappa_sca_m = np.zeros(len(lam))
    for i in np.arange(len(lam)):
        kappa_abs_m[i] = sum(np.transpose(kappa_abs[:, i]) * sig, 0)
        kappa_sca_m[i] = sum(np.transpose(kappa_sca[:, i]) * sig, 0)
    #
    # ... and return them
    #
    return kappa_abs_m, kappa_sca_m


def get_mie_coefficients(A, LAM, diel_constants, bhmie_function=bhmie_function, nang=3, return_all=False):
    """
    This calculates the opacity for the given dielectric constants for all
    grain sizes and wavelength specified in LAM and A.

    Arguments:
    ----------
    A : array
    : all the grain sizes for which the opacities are calculated

    LAM : array
    :    all the wavelength in cm at which the opacities are calculated

    diel_constants : object of class diel_const
    :    the dielectric constants that are used for the calculation

    Keywords:
    ---------

    method: callable
    : a function that carries out the Mie calculation with this signature
        S1, S2, Qext, Qabs, Qsca, Qback, gsca = bhmie_function(x, (n, k), n_angles)

    nang: int
    :    Number of angles for scattering function (not supported by all methods)

    return_all : bool
    :    wether or not to include also the asymmetry factor or scattering function in the output.
         If true, the output elements 2 and onwards are gg_sca, s_1, s_2

    Output:
    -------
    q_abs,q_sca : array

    if return_all is true, the output is
    q_abs,q_sca, gg_sca, s_1, s_2
    """
    from numpy import zeros, maximum
    from scipy.optimize import fsolve
    #
    # feed the bhmie function
    # use the first entries
    #
    q_abs = zeros([len(A), len(LAM)])
    q_sca = zeros([len(A), len(LAM)])
    gg_sca = zeros([len(A), len(LAM)])
    s_1 = zeros([len(A), len(LAM), nang], dtype=complex)
    s_2 = zeros([len(A), len(LAM), nang], dtype=complex)
    NMXX = 200000
    #
    # the wave length loop
    #
    for ilam, lam in enumerate(LAM):
        #
        # display progress
        #
        progress_bar((ilam + 1.0) / len(LAM) * 100, 'Mie')
        #
        # interpolate the refr. index
        #
        n, k = diel_constants.nk(lam)
        #
        # define the size parameter
        #
        X = 2. * np.pi / lam * A
        #
        # define the cutoff where no convergence is reached
        #
        Y = X * (n + k * 1j)
        Y = abs(Y)
        Xstop = X + 4. * X**.333333 + 2.0
        nmx = maximum(Xstop, Y).astype(int) + 15
        mask = nmx < NMXX
        #
        # cut it such that only converging terms are included
        # and extrapolate the missing parts (see below)
        #
        X_cut = X[mask]
        if len(X_cut) == 0:
            def f(x):
                return x + 4. * x**0.33333 + 2.0 - NMXX
            x_max = fsolve(f, NMXX)  # /30.
            X_cut = [x_max]
        #
        # loop through the sizes
        #
        for ia, x in enumerate(X_cut):
            S1, S2, Qext, Qabs, Qsca, Qback, gsca = bhmie_function(x, complex(n, k), nang)
            q_abs[ia, ilam] = Qabs
            q_sca[ia, ilam] = Qsca
            gg_sca[ia, ilam] = gsca.real
            s_1[ia, ilam, :] = S1[:nang]
            s_2[ia, ilam, :] = S2[:nang]
        #
        # extrapolate for large grains
        #
        q_abs[ia + 1:, ilam] = q_abs[ia, ilam]
        q_sca[ia + 1:, ilam] = q_sca[ia, ilam]
        gg_sca[ia + 1:, ilam] = gg_sca[ia, ilam]

    if return_all:
        return q_abs, q_sca, gg_sca, s_1, s_2
    else:
        return q_abs, q_sca


def compare_nk(constants, lmin=1e-5, lmax=1e3, orig_data=False, ax=None):
    """
    Compares the dielectric functions c1 and c2 (their n and k values)
    by plotting them on the range from amin to amax.

    Keywords:
    ---------

    constants : list of instances of diel_constant
        the dielectric constants to compare

    lmin : float
    :    lower bound of the plotting range

    lmax : float
    :    upper bound of the plotting range

    Keywords:
    ---------

    orig_data : bool
        if true, then just plot the original data of each object

    ax : plt.axes
        plot into those axes, or create new figure/axes

    Output:
    -------
    A plot comparing both optical property functions
    """
    import matplotlib.pyplot as plt

    nlam = 100
    lam = np.logspace(np.log10(lmin), np.log10(lmax), nlam)

    if ax is None:
        f, ax = plt.subplots()

    for n, c in enumerate(constants):
        nk = np.zeros([len(lam), 2])

        if orig_data and (c._n is not None) and (c._l is not None):
            lam = c._l
            nk = np.array([c._n, c._k]).T
        else:
            for i in np.arange(nlam):
                try:
                    nk[i] = c.nk(lam[i])
                except BaseException:
                    nk[i] = np.nan

        mask = np.invert(np.isnan(nk[:, 0]))
        ax.loglog(lam[mask], nk[:, 0][mask], label='$n_1$ - {}'.format(c.material_str), c=f'C{n}', alpha=0.7, ls='-')
        mask = np.invert(np.isnan(nk[:, 1]))
        ax.loglog(lam[mask], nk[:, 1][mask], label='$k_1$ - {}'.format(c.material_str), c=f'C{n}', alpha=0.7, ls='--')

    ax.set_xlabel('wavelength [cm]')
    ax.set_ylabel('$n, k$')
    ax.legend(loc='best')

    return ax


def get_default_diel_constants(extrapol=False, lmax=None):
    """
    This method calculates the mixed mie coefficients as in Ricci et al. 2010.

    The densities and volume fractions stated in Ricci et al. 2010 contained
    typos. The values used here come from L. Ricci, private communications.

    |                 | silicates |  carbonaceous | water ice | vacuum |
    |:---------------:|:---------:|:-------------:|:---------:|:------:|
    | volume fraction | 0.07      |     0.21      |    0.42   | 0.30   |
    | solid densities | 3.5 g/cc  |    2.5 g/cc   |  1 g/cc   | 0 g/cc |

    Keywords:
    ---------

    extrapol : bool
        whether the optical constants do extrapolation or not

    Output:
    -------

    diel_constants : object of class diel_const
        the mixed dielectric constants

    rho_s : float
        the material density of the particles in g/cm**3
    """
    if extrapol and (lmax is None):
        raise ValueError('need to set lmax if extrapol is True')

    c1 = diel_draine2003_astrosil()
    c2 = diel_zubko_carbon(extrapol=extrapol, lmax=lmax)
    c3 = diel_warren(new=True)
    c4 = diel_vacuum()

    # c1 = diel_luca('silicate', extrapol=extrapol, lmin=1e-6, lmax=100.0)
    # c2 = diel_luca('carbon', extrapol=extrapol, lmin=1e-6, lmax=100.0)
    # c3 = diel_luca('ice', extrapol=extrapol, lmin=1e-6, lmax=100.0)
    # c4 = diel_vacuum()
    constants = [c1, c2, c3, c4]

    # after Lucas thesis, the fractions in Ricci+2010 are typos

    vol_fract = [0.07, 0.21, 0.42, 0.30]
    densities = [3.50, 2.50, 1.00, 0.00]

    diel_constants = diel_mixed(constants, vol_fract, rule='Bruggeman')
    rho_s = sum(densities * np.array(vol_fract))

    return diel_constants, rho_s


def get_opacities(a, lam, rho_s=None, diel_const=None, bhmie_function=bhmie_function,
                  return_all=False, extrapol=False):
    """
    Calculates opacities according to some specified method for
    a given size- and wavelength grid. If diel_const and rho_s is not given,
    default opacities are used. Otherwise both rho_s and diel_const are needed.

    Arguments:
    ----------

    a : array
        The grain size grid in cm

    lam : array
        the wavelength grid in cm

    Keywords:
    ---------

    diel_const : dielectric constant
        in case other dielectric constants should be used

    rho_s : float
        material density of each size [g/cm^3]

    bhmie_function : callable
        which function to use for the mie calculation

    return_all : bool
        default False: return just kappa_abs, kappa_sca, rho_s
        True: return kappa_*, assymetry factor, S1, S2, and rho_s

    extrapol : bool
        whether to extrapolate default optical constants if lam is outside the
        wavelength range of the data

    Output:
    -------
    kappa_abs, kappa_sca : arrays
        absorption and scattering opacities [g/cm^2]

    gg_sca : array

    s_1, s_2 : arrays
        the Mueller matrix elements

    rho_s : float
        material density of the grains [g/cm^3]

    """
    if (diel_const is None and rho_s is not None) or (diel_const is not None and rho_s is None):
        raise AssertionError('diel_const and rho_s need to be both given or both be None')

    if diel_const is None and rho_s is None:
        if extrapol:
            print('Note: wavelength range outside data. Extrapolation will be done but is uncertain.')
        diel_const, rho_s = get_default_diel_constants(extrapol=extrapol, lmax=lam[-1])

    m = 4 * np.pi / 3. * rho_s * a**3

    q_abs, q_sca, gg_sca, s_1, s_2 = get_mie_coefficients(a, lam, diel_const, return_all=True, bhmie_function=bhmie_function)

    kappa_abs, kappa_sca = get_kappa_from_q(a, m, q_abs, q_sca)

    if return_all:
        return kappa_abs, kappa_sca, gg_sca, s_1, s_2, rho_s
    else:
        return kappa_abs, kappa_sca, rho_s
