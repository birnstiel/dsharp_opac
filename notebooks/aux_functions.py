from astropy import units as u
from astropy import constants as c
import matplotlib.pyplot as plt
import numpy as np

# style settings for DSHARP

style = [
    'default',
    {
        'figure.figsize': (3.5, 3.5 / 1.618),
        'font.size': 9,  # 12
        'image.cmap': 'inferno',
        'figure.dpi': 200,
        'font.family': ['Times', 'Times New Roman'],
        'xtick.top': True,
        'xtick.direction': 'in',
        'ytick.right': True,
        'ytick.direction': 'in',
        'mathtext.fontset': 'cm'
        }]


def set_style():
    "Set the style of DSHARP"
    plt.style.use(style)


def get_style():
    "Get a style context manager for DSHARP"
    return plt.style.context(style)


def print_ratio(f, r0=1.618):
    """
    Call this function after each plot to check the aspect ratio.
    """
    for ax in f.get_axes():
        r = ax.spines['bottom'].get_extents().width / ax.spines['left'].get_extents().height
        print('r = {}; f = {}'.format(r, r0 / r))


def planck_B_nu(freq, T):
    """"
    Calculates the value of the Planck-Spectrum
    B(nu,T) of a given frequency nu and temperature T

    Arguments
    ---------
    nu : float or array
        frequency in 1/s or with astropy.units

    T: float
        temperature in K or in astropy.units

    Returns:
    --------
    B : float
        value of the Planck-Spectrum at frequency nu and temperature T
        units are using astropy.units if the input values use those, otherwise
        cgs units: erg/(s*sr*cm**2*Hz)

    """
    if isinstance(T, u.quantity.Quantity):
        use_units = True
    else:
        T = T * u.K
        use_units = False

    if not isinstance(freq, u.quantity.Quantity):
        freq *= u.Hz

    T = np.array(T.value, ndmin=1) * T.unit
    freq = np.array(freq.value, ndmin=1) * freq.unit

    f_ov_T = freq[np.newaxis, :] / T[:, np.newaxis]
    mx = np.floor(np.log(np.finfo(f_ov_T.ravel()[0].value).max))
    exp = np.minimum(f_ov_T * c.h / c.k_B, mx)
    exp = np.maximum(exp, -mx)

    output = 2 * c.h * freq**3 / c.c**2 / (np.exp(exp) - 1.0) / u.sr

    cgsunit = 'erg/(s*sr*cm**2*Hz)'
    if use_units:
        return output.to(cgsunit).squeeze()
    else:
        return output.to(cgsunit).value.squeeze()


def planck_dBnu_dT(freq, T):
    """This function computes the temperature derivative of the blackbody function

            dB_nu(T)        2 h^2 nu^4      exp(h nu / kT)        1
            --------   =    ---------- ------------------------  ---
               dT            k c^2    [ exp(h nu / kT) - 1 ]^2  T^2

    Arguments
    ---------

    freq : float or array
        frequency [in Hz or with astropy.units]

    temp : float or array
        temperature [in K or with astropy.units]

    Returns:
    --------
    dBdT: the temperature derivative of the planck spectrum
        units are using astropy.units if the input values use those, otherwise
        cgs units: erg/(K*s*sr*cm**2*Hz)

    """
    if isinstance(T, u.quantity.Quantity):
        use_units = True
    else:
        T = T * u.K
        use_units = False

    if not isinstance(freq, u.quantity.Quantity):
        freq = freq * u.Hz

    T = np.array(T.value, ndmin=1) * T.unit
    freq = np.array(freq.value, ndmin=1) * freq.unit

    f_ov_T = freq[np.newaxis, :] / T[:, np.newaxis]
    mx = np.floor(np.log(np.finfo(f_ov_T.ravel()[0].value).max))
    exp = np.minimum(f_ov_T * c.h / c.k_B, mx)
    exp = np.maximum(exp, -mx)
    exp = np.exp(exp)

    dBdT = 2 * c.h**2 * freq[np.newaxis, :]**4 / (c.k_B * c.c**2 * T[:, np.newaxis]**2) * \
        1. / (exp * (1. - 1. / exp)**2) / u.sr

    cgsunit = 'erg/(s*Hz*sr*cm**2*K)'
    if use_units:
        return dBdT.to(cgsunit).squeeze()
    else:
        return dBdT.to(cgsunit).value.squeeze()


def permutate_fiducial(params):
    """
    This function creates all permutations of parameters where only one parameter
    is allowed to change at a time. The fiducial model is the central entry of
    each array, or the one left of the center if there is an even number of
    entries. The output will be an array with (number of input arrays) columns.
    At least 2 arrays need to be given.

    Arguments:
    ----------
    a: list
    :    list of parameter lists, each containing the allowed parameter values

    Example
    >>> permutate_fiducial([[1],[1,2],[1,2,3]])
    array([[ 1.,  2.,  2.],
           [ 1.,  1.,  1.],
           [ 1.,  1.,  3.],
           [ 1.,  1.,  2.]])
    """
    def get_fiducial_idx(a):
        return int(np.ceil(len(a) / 2.) - 1)
    param_sets = np.zeros([sum([len(i) - 1 for i in params]) + 1, len(params)])
    for i, param in enumerate(params):
        param_sets[:, i] = param[get_fiducial_idx(param)]
    k = 0
    for i, param in enumerate(params):
        for j, value in enumerate(param):
            if j == get_fiducial_idx(param):
                continue
            else:
                param_sets[k, i] = value
                k += 1
    return param_sets


def t_sat_water(Sigma, M_star, r, f_h2o=0.005):
    """Calculate water sublimation temperature using values from Leger et al. 1985.

    Parameters
    ----------
    Sigma : float
        total gas surface density [g/cm^2]

    M_star : float
        stellar mass [g]

    r : float
        radius [au]

    f_h2o : float
        water abundance

    Returns
    -------
    float
        sublimation temperature

    """
    from scipy.optimize import fsolve

    torr = 101325. / 760. * (1. * u.Pa).cgs.value
    muw  = 18.01528
    mug  = 2.3
    k_b  = c.k_B.cgs.value
    m_p  = c.m_p.cgs.value
    p_0  = 1.9e10 * torr
    dH   = 6070.
    om   = np.sqrt(c.G.cgs.value * M_star / r**3)

    A = f_h2o * Sigma * om / p_0 * np.sqrt(k_b * mug / (2 * np.pi * muw**2 * m_p))

    def fct(T):
        return A * np.sqrt(T) - np.exp(-dH / T)

    return fsolve(fct, 170)[0]
