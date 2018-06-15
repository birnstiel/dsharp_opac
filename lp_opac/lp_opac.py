#!/usr/bin/env python
"""
This module contains opacity scripts and all the helper and testing routines.

- dielectric functions are objects, see bhmie.diel_*
  the functions `diel_*.nk` return the optical properties
"""
import numpy as np
import os
import sys
import pkg_resources

from .bhmie import bhmie

pi = np.pi

try:
    import pymiecoated

    def bhmie_pymiecoated(x, nk, angles):
        """
        Wrapper to the Mie code of `pymiecoated`

        x : float
            size parameter x = 2 pi a / lambda

        nk : complex
            complex refractive index n + i k

        angles : array
            array of angles in degree for which to return S1 & S2

        Output:
        -------
        S1, S2, Qext, Qsca, Qback, gsca

        S1, S2 : arrays
            the matrix elements as function of angle

        Qext, Qsca, Qback : float
            the extinction, scattering backscatterin coefficients

        gsca : float
            Henyey-Greenstein asymmetry factor
        """
        #
        # the other code
        #
        opac = pymiecoated.Mie(x=x, m=nk)
        Qext = opac.qext()
        Qsca = opac.qsca()
        Qback = opac.qratio()
        gsca = opac.asy()

        S1 = np.zeros(len(angles))
        S2 = np.zeros(len(angles))

        for i, angle in enumerate(angles):
            S1[i], S2[i] = opac.S12(np.cos(angle / 180 * np.pi))

        return S1, S2, Qext, Qsca, Qback, gsca
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
            raise NameError('{}: wavelength {:g} outside data-range [{:g},{g}]'.format(type(self).__name__, l, self._lmin, self._lmax))
        return 10.**np.interp(np.log10(l), self._ll, self._ln), 10.**np.interp(np.log10(l), self._ll, self._lk)


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


class diel_wd03_sil(diel_const):
    """
    Returns the dielectric constants for astronomical silicates from
    Weingartner & Draine 2003. The data comes from callindex.out_sil.D03`
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
        self.material_str = 'Astronomical Silicates (Weingartner & Draine 2003)'
        self.datafile = pkg_resources.resource_filename(__name__, os.path.join(
            'optical_constants', 'draine', 'callindex.out_silD03'))
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


class diel_dl84_astrosil(diel_const):
    """
    Returns the dielectric constants for astronomical silicates from
    Draine & Lee 1984. The data comes from eps_Sil
    which was downloaded from

        ftp://ftp.astro.princeton.edu/draine/dust/diel/eps_Sil

    on 2014-10-28--15:24 EDT
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
        self.headerinfo = [f.readline() for i in range(6)]
        data = np.loadtxt(f)
        #
        # assign wavelength and optical constants
        #
        self._l = data[::-1, 0] * 1e-4
        self._n = data[::-1, 1] + 1.
        self._k = data[::-1, 2]
        self._ll = np.log10(self._l)
        self._ln = np.log10(self._n)
        self._lk = np.log10(self._k)
        self._lmin = self._l.min()
        self._lmax = self._l.max()


class diel_vacuum(diel_const):
    """
    Returns the dielectric constants for vacuum
    """

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
    al. 1996 (the BE values). The data was OCRed by Til Birnstiel, no
    guarantee for correctness.

    Keywords:
    ---------

    extrapol : bool
    :    whether or not to extrapolate beyond ~0.19 mm

    lmax : float
    :    upper limit for extrapolation range
    """

    def __init__(self, extrapol=False, lmax=1.0):
        """
        Overwrite the initialization of the parent class
        """
        #
        # open files, read data
        #
        self.material_str = 'Carbonaceous Grains (Zubko et al. 1996)'
        self.datafile = pkg_resources.resource_filename(
            __name__, os.path.join('optical_constants', 'zubko+1996'))

        E = np.loadtxt(pkg_resources.resource_filename(
            __name__, os.path.join(self.datafile, 'zubko_E.txt')))[-1::-1]
        n = np.loadtxt(pkg_resources.resource_filename(
            __name__, os.path.join(self.datafile, 'zubko_n_BE.txt')))[-1::-1]
        k = np.loadtxt(pkg_resources.resource_filename(
            __name__, os.path.join(self.datafile, 'zubko_k_BE.txt')))[-1::-1]
        l = 0.00012398419292004205 / E  # E = h*c/lambda in CGS # noqa
        #
        # extrapolate
        #
        if extrapol:
            le = np.logspace(np.log10(l[-1]), np.log10(lmax), 10)
            from scipy.optimize import curve_fit

            def f(x, a, b, c):
                return a + b * x + c * x**2

            i_min = abs(l - 1e-2).argmin()
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


class diel_warren(diel_const):
    """
    Returns the dielectric constants for water ice according
    to the data by Warren 1986 or the newer version

    the new, updated data was downloaded from

        http://www.atmos.washington.edu/ice_optical_constants/

    on Jan 23, 2014.

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
        self.material_str = 'Water Ice (Warren, %s data)' % (new * 'new' + (not new) * 'old')
        #
        # set the file name
        #
        if new:
            fname = 'IOP_2008_ASCIItable.dat'
        else:
            fname = 'warren_1986.txt'
        self.datafile = pkg_resources.resource_filename(
            __name__, os.path.join('optical_constants', 'warren', fname))
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


def get_total_opacity(a, lam, n, rho_s, diel_constants, q_abs=None, q_sca=None):
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

    diel_constants : object of class diel_const
        the dielectric constants to be used

    Keywords:
    ---------

    q_abs,q_sca : array
        if the opacity for all sizes and wavelength has already
        been calculated, then you can pass it along, otherwise
        it will be calculated on the fly.

    Output:
    -------
    kappa_abs,kappa_sca : array
        the opacity at each wavelength averaged over the size
        distribution and normalized be per 1 g of dust.
    """
    from numpy import zeros, log
    #
    # some basic conversions
    #
    m = 4. * pi / 3. * rho_s * a**3
    C = log(a[2] / a[1])
    sig = n * m * a * C
    sig = sig / sum(sig)
    #
    # calculate the opacities ...
    #
    if q_abs is None or q_sca is None:
        q_abs, q_sca = get_mie_coefficients(a, lam, diel_constants)
    kappa_abs = q_abs * np.tile(pi * a**2 / m, [len(lam), 1]).transpose()
    kappa_sca = q_sca * np.tile(pi * a**2 / m, [len(lam), 1]).transpose()
    #
    # ... average them over the size distribution ...
    #
    kappa_abs_m = zeros(len(lam))
    kappa_sca_m = zeros(len(lam))
    for i in np.arange(len(lam)):
        kappa_abs_m[i] = sum(np.transpose(kappa_abs[:, i]) * sig, 0)
        kappa_sca_m[i] = sum(np.transpose(kappa_sca[:, i]) * sig, 0)
    #
    # ... and return them
    #
    return kappa_abs_m, kappa_sca_m


def get_opacity_from_distribution(a, lam, n, dc, rho_s, bhmie_function=bhmie):
    """
    A simple wrapper for the opacity functions: turns a given
    dielectric constant and a given size distribution into the
    size-averaged opacities on a given wavelength grid

    Arguments:
    ----------
    a : array
    :    grain size grid in cm

    lam : array
    :    wavelength grid in cm

    n : array
    :    grain size distribution n(a) in cm^-3

    dc : instance of class diel_const
    :    the n and k values to use

    rho_s : float
    :    the material density of the grains

    Keywords:
    ---------
    method : str
    :     which of the bhmie codes to use

    Returns:
    --------
    kap_abs_t,kap_sca_t : arrays
    :    the absorption and scattering opacity in cm^2/g

    """
    q_abs, q_sca = get_mie_coefficients(a, lam, dc, bhmie_function=bhmie_function)
    kap_abs_t, kap_sca_t = get_total_opacity(a, lam, n, rho_s, None, q_abs=q_abs, q_sca=q_sca)
    return kap_abs_t, kap_sca_t


def get_mie_coefficients(A, LAM, diel_constants, bhmie_function=bhmie, nang=3, return_all=False):
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
        S1, S2, Qext, Qsca, Qback, gsca = bhmie(x, complex(n, k), angles)

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
    s_1 = zeros([len(A), len(LAM), 2 * nang - 1], dtype=complex)
    s_2 = zeros([len(A), len(LAM), 2 * nang - 1], dtype=complex)
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
        X = 2. * pi / lam * A
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
            S1, S2, Qext, Qsca, Qback, gsca = bhmie_function(x, complex(n, k), np.linspace(0., 180., nang))
            q_abs[ia, ilam] = Qext - Qsca
            q_sca[ia, ilam] = Qsca
            gg_sca[ia, ilam] = gsca
            s_1[ia, ilam, :] = S1
            s_2[ia, ilam, :] = S2
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


def compare_nk(c1=None, c2=None, amin=1e-5, amax=1e3):
    """
    Compares the dielectric functions c1 and c2 (the n and k values)
    by plotting them on the range from amin to amax.

    Keywords:
    ---------

    c1 : instance of diel_constant
    :    first diel. function, defaults to ice without extrapolation

    c2 : instance of diel_constant
    :    second diel. function, defaults to ice with extrapolation

    amin : float
    :    lower bound of the plotting range

    amax : float
    :    upper bound of the plotting range

    Output:
    -------
    A plot comparing both optical property functions
    """
    import matplotlib.pyplot as plt

    if c1 is None:
        c1 = diel_luca('ice', extrapol=False)
    if c2 is None:
        c2 = diel_luca('ice', extrapol=True, lmax=amax)

    na = 100

    a = np.logspace(np.log10(amin), np.log10(amax), na)

    nk1 = np.zeros([len(a), 2])
    nk2 = np.zeros([len(a), 2])

    for i in np.arange(na):
        try:
            nk1[i] = c1.nk(a[i])
        except BaseException:
            pass
        try:
            nk2[i] = c2.nk(a[i])
        except BaseException:
            pass

    f, ax = plt.subplots()
    ax.loglog(a, nk1[:, 0], label='$n_1$%s' % c1.material_str, c='C0', alpha=0.7, ls='-')
    ax.loglog(a, nk2[:, 0], label='$n_2$%s' % c1.material_str, c='C0', alpha=0.7, ls='--')

    if c1.extrapol is True:
        ax.axvline(c1._lmin, c='C0')
        ax.axvline(c1._lmax, c='C0')

    ax.loglog(a, nk1[:, 1], label='$k_1$%s' % c1.material_str, c='C1', alpha=0.7, ls='-')
    ax.loglog(a, nk2[:, 1], label='$k_2$%s' % c1.material_str, c='C1', alpha=0.7, ls='--')

    if c2.extrapol is True:
        ax.axvline(c2._lmin, c='C1')
        ax.axvline(c2._lmax, c='C1')

    ax.set_xlabel('wavelength')
    ax.set_ylabel('$n, k$')
    ax.legend(loc='best')


def get_default_diel_constants():
    """
    This method calculates the mixed mie coefficients as in Ricci et al. 2010.

    The densities and volume fractions stated in Ricci et al. 2010 contained
    typos. The values used here come from L. Ricci, private communications.

    |                 | silicates |  carbonaceous | water ice | vacuum |
    |:---------------:|:---------:|:-------------:|:---------:|:------:|
    | volume fraction | 0.07      |     0.21      |    0.42   | 0.30   |
    | solid densities | 3.5 g/cc  |    2.5 g/cc   |  1 g/cc   | 0 g/cc |

    Output:
    -------

    diel_constants : object of class diel_const
        the mixed dielectric constants

    rho_s : float
        the material density of the particles in g/cm**3
    """
    # c1 = diel_wd03_sil()
    # c2 = diel_zubko_carbon(extrapol=True)
    # c3 = diel_warren(new=False)
    # c4 = diel_vacuum()

    c1 = diel_luca('silicate', extrapol=True, lmax=100.0)
    c2 = diel_luca('carbon', extrapol=True, lmax=100.0)
    c3 = diel_luca('ice', extrapol=True, lmax=100.0)
    c4 = diel_vacuum()
    constants = [c1, c2, c3, c4]

    # after Lucas thesis, the fractions in Ricci+2010 are typos

    vol_fract = [0.07, 0.21, 0.42, 0.30]
    densities = [3.50, 2.50, 1.00, 0.00]

    diel_constants = diel_mixed(constants, vol_fract, rule='Bruggeman')
    rho_s = sum(densities * np.array(vol_fract))

    return diel_constants, rho_s


def get_default_opacities(a, lam, bhmie_function=bhmie, return_all=False):
    """
    Calculates opacities according to some specified method for
    a given size- and wavelength grid.

    Arguments:
    ----------

    a : array
        The grain size grid in cm

    lam : array
        the wavelength grid in cm

    Keywords:
    ---------

    bhmie_function : callable
        which function to use for the mie calculation

    return_all : bool
        default False: return just kappa_abs, kappa_sca, rho_s
        True: return kappa_*, assymetry factor, S1, S2, and rho_s

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
    diel_const, rho_s = get_default_diel_constants()

    m = 4 * np.pi / 3. * rho_s * a**3

    q_abs, q_sca, gg_sca, s_1, s_2 = get_mie_coefficients(a, lam, diel_const, return_all=True, bhmie_function=bhmie_function)

    kappa_abs = q_abs * np.tile(np.pi * a**2 / m, [len(lam), 1]).transpose()
    kappa_sca = q_sca * np.tile(np.pi * a**2 / m, [len(lam), 1]).transpose()

    if return_all:
        return kappa_abs, kappa_sca, gg_sca, s_1, s_2, rho_s
    else:
        return kappa_abs, kappa_sca, rho_s
