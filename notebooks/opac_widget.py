import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.widgets as w

import astropy.constants as c

import aux_functions as aux
import dsharp_opac as opacity

try:
    import coala
    coala_installed = True
except ImportError:
    coala_installed = False

aux.set_style()

# ------------------------
# GENERAL CALCULATION PART
# ------------------------

# load default opacities

with np.load(opacity.get_datafile('default_opacities_smooth.npz')) as d:
    a_w     = d['a']
    gsca_w  = d['g']
    lam_w   = d['lam']
    k_abs_w = d['k_abs']
    k_sca_w = d['k_sca']
    k_ext_w = k_abs_w + (1 - gsca_w) * k_sca_w

# load dry opacities

with np.load(opacity.get_datafile('icefree_opacities_smooth.npz')) as d:
    a_d      = d['a']
    gsca_d   = d['g']
    lam_d    = d['lam']
    k_abs_d  = d['k_abs']
    k_sca_d  = d['k_sca']
    k_ext_d = k_abs_d + (1 - gsca_d) * k_sca_d

# sanity check

if np.allclose(lam_w, lam_d) and np.allclose(a_w, a_d):
    lam = lam_w
    a = a_w
else:
    raise ValueError('grids of the opacities do not match')

# get constants and calculate derived quantities

nu    = c.c.cgs.value / lam
au    = c.au.cgs.value
M_sun = c.M_sun.cgs.value

# set initial parameter values

r       = au
sigma_g = 200.0
d2g     = 0.01
rho_s   = 1.6686
M_star  = M_sun
v_frag  = 100.0
alpha   = 0.001
T       = 100
mu      = 2.3

# set parameter names and ranges

par_names = ['$T$', r'$v_\mathrm{frag}$', r'$\alpha$', r'$\Sigma_\mathrm{g}$']
par_start = [1,    50,     1e-5,  1]        # noqa - minimum parameter value
par_end   = [1500, 3000,   1e-1,  1e4]      # noqa - maximum parameter value
par_ini   = [T,    v_frag, alpha, sigma_g]  # noqa - initial parameter value
n_grid    = [10,   10,     10,    10]       # noqa - number of parameter values
n_params  = len(par_names)
par_grid = [
    np.logspace(np.log10(p_start), np.log10(p_end), _n)
    for p_start, p_end, _n in zip(par_start, par_end, n_grid)
    ]

# define the power-law

q = 3.5
amax = 0.1
power_law = a**(4 - q)
power_law[a > amax] = 0
power_law = power_law / power_law.sum()

# assuming T is the first parameter: calculate
# the quantities for the power-law distribution

k_P_pf   = np.zeros([n_grid[0]])
k_R_pf   = np.zeros([n_grid[0]])

k_abs_p_w = (k_abs_w.T * power_law[None, :]).sum(1)
k_ext_p_w = (k_ext_w.T * power_law[None, :]).sum(1)

# the dry ones

k_abs_p_d = (k_abs_d.T * power_law[None, :]).sum(1)
k_ext_p_d = (k_ext_d.T * power_law[None, :]).sum(1)

for it, _T in enumerate(par_grid[0]):
    Bnu    = aux.planck_B_nu(nu, _T)
    dBnudT = aux.planck_dBnu_dT(nu, _T)
    B      = np.trapz(Bnu, x=nu)
    dBdT   = np.trapz(dBnudT, x=nu)

    if _T < 170:
        k_P_pf[it]  = np.trapz(Bnu * k_abs_p_w, x=nu) / B
        k_R_pf[it]  = dBdT / np.trapz(dBnudT / k_ext_p_w, x=nu)
    else:
        k_P_pf[it]  = np.trapz(Bnu * k_abs_p_d, x=nu) / B
        k_R_pf[it]  = dBdT / np.trapz(dBnudT / k_ext_p_d, x=nu)

# --------------------
# DEFINE DATA FUNCTION
# --------------------


def get_data(values):

    T, v_frag, alpha, sigma_g = values

    # get the size distribution fit number 1

    f1 = opacity.get_B11_fit(
        T, a, r=r, sigma_g=sigma_g, d2g=d2g, rho_s=rho_s, M_star=M_star,
        v_frag=v_frag, alpha=alpha)
    f1 = f1 / f1.sum()

    # get the size distribution fit number 2

    f2, amax = opacity.get_B11S_fit(
        T, a, r=r, sigma_g=sigma_g, d2g=d2g, rho_s=rho_s, M_star=M_star,
        v_frag=v_frag, alpha=alpha)
    f2 = (f2 * a) / (f2 * a).sum()

    if T < aux.t_sat_water(sigma_g, M_star, r):
        # sum the absorption opacity

        k_abs_f1 = (k_abs_w.T * f1[None, :]).sum(1)
        k_abs_f2 = (k_abs_w.T * f2[None, :]).sum(1)
        k_abs_pl = (k_abs_w.T * power_law[None, :]).sum(1)

        # sum the EXTINCTION opacity

        k_ext_f1 = (k_ext_w.T * f1[None, :]).sum(1)
        k_ext_f2 = (k_ext_w.T * f2[None, :]).sum(1)
        k_ext_pl = (k_ext_w.T * power_law[None, :]).sum(1)
    else:
        # sum the DRY absorption opacity

        k_abs_f1 = (k_abs_d.T * f1[None, :]).sum(1)
        k_abs_f2 = (k_abs_d.T * f2[None, :]).sum(1)
        k_abs_pl = (k_abs_d.T * power_law[None, :]).sum(1)

        # sum the DRY EXTINCTION opacity

        k_ext_f1 = (k_ext_d.T * f1[None, :]).sum(1)
        k_ext_f2 = (k_ext_d.T * f2[None, :]).sum(1)
        k_ext_pl = (k_ext_d.T * power_law[None, :]).sum(1)

    # calculate Planck opacity for the fit and the power-law

    Bnu    = aux.planck_B_nu(nu, T)
    dBnudT = aux.planck_dBnu_dT(nu, T)
    B      = np.trapz(Bnu, x=nu)

    k_P_f1 = np.trapz(Bnu * k_abs_f1, x=nu) / B
    k_P_f2 = np.trapz(Bnu * k_abs_f2, x=nu) / B
    k_P_pl = np.trapz(Bnu * k_abs_pl, x=nu) / B

    # calculate Rosseland opacity for the fit and the power-law

    dBdT   = np.trapz(dBnudT, x=nu)
    k_R_f1 = dBdT / np.trapz(dBnudT / k_ext_f1, x=nu)
    k_R_f2 = dBdT / np.trapz(dBnudT / k_ext_f2, x=nu)
    k_R_pl = dBdT / np.trapz(dBnudT / k_ext_pl, x=nu)

    return {
        'f1': f1,
        'f2': f2,
        'amax': amax,
        'Bnu': Bnu,
        'k_abs_f1': k_abs_f1,
        'k_abs_f2': k_abs_f2,
        'k_P_f1': k_P_f1,
        'k_P_f2': k_P_f2,
        'k_P_pl': k_P_pl,
        'k_R_f1': k_R_f1,
        'k_R_f2': k_R_f2,
        'k_R_pl': k_R_pl,
        }


# This part calculates the k_P for all combinations
# this takes quite a while if n_grid > 10, so we
# store it and load it if if it's available

fname = 'kPR2.npz'
if os.path.isfile(fname):
    k_P_f2_array = np.load(fname)['k_P_f2_array']
    k_R_f2_array = np.load(fname)['k_R_f2_array']
    k_P_pl_array = np.load(fname)['k_P_pl_array']
    k_R_pl_array = np.load(fname)['k_R_pl_array']
    amax_array   = np.load(fname)['amax_array']
    T_array      = np.load(fname)['T_array']
else:
    k_P_f2_array = np.zeros([len(p) for p in par_grid])
    k_R_f2_array = np.zeros([len(p) for p in par_grid])
    k_P_pl_array = np.zeros([len(p) for p in par_grid])
    k_R_pl_array = np.zeros([len(p) for p in par_grid])
    amax_array   = np.zeros([len(p) for p in par_grid])

    for i1 in range(n_grid[0]):
        for i2 in range(n_grid[1]):
            for i3 in range(n_grid[2]):
                for i4 in range(n_grid[3]):
                    res = get_data([
                        par_grid[0][i1],
                        par_grid[1][i2],
                        par_grid[2][i3],
                        par_grid[3][i4],
                        ])

                    k_R_f2_array[i1, i2, i3, i4] = res['k_R_f2']
                    k_P_f2_array[i1, i2, i3, i4] = res['k_P_f2']
                    k_R_pl_array[i1, i2, i3, i4] = res['k_R_pl']
                    k_P_pl_array[i1, i2, i3, i4] = res['k_P_pl']
                    amax_array[i1, i2, i3, i4]   = res['amax']

    T_array = par_grid[0]
    np.savez_compressed(fname, **{
        'k_P_f2_array': k_P_f2_array,
        'k_R_f2_array': k_R_f2_array,
        'k_P_pl_array': k_P_pl_array,
        'k_R_pl_array': k_R_pl_array,
        'amax_array': amax_array,
        'T_array': T_array})

# find the max & min values

kP2min  = k_P_f2_array.copy()
kP2max  = k_P_f2_array.copy()

kR2min  = k_R_f2_array.copy()
kR2max  = k_R_f2_array.copy()

kPplmin  = k_P_pl_array.copy()
kPplmax  = k_P_pl_array.copy()

kRplmin  = k_R_pl_array.copy()
kRplmax  = k_R_pl_array.copy()

kP2min[np.isnan(kP2min)] = k_P_f2_array.max()
kP2max[np.isnan(kP2max)] = k_P_f2_array.min()
kR2min[np.isnan(kR2min)] = k_R_f2_array.max()
kR2max[np.isnan(kR2max)] = k_R_f2_array.min()

kPplmin[np.isnan(kPplmin)] = k_P_pl_array.max()
kPplmax[np.isnan(kPplmax)] = k_P_pl_array.min()
kRplmin[np.isnan(kRplmin)] = k_R_pl_array.max()
kRplmax[np.isnan(kRplmax)] = k_R_pl_array.min()

for i in range(n_params - 1):
    kP2min = np.nanmin(kP2min, -1)
    kP2max = np.nanmax(kP2max, -1)
    kR2min = np.nanmin(kR2min, -1)
    kR2max = np.nanmax(kR2max, -1)

    kPplmin = np.nanmin(kPplmin, -1)
    kPplmax = np.nanmax(kPplmax, -1)
    kRplmin = np.nanmin(kRplmin, -1)
    kRplmax = np.nanmax(kRplmax, -1)

# -----------------
# INITIALIZE FIGURE
# -----------------

# set plot spacing

hslider  = 0.03
dyslider = hslider + 0.01
xslider  = 0.3
wslider  = 0.3
panelbot = 0.0
controlh = panelbot + n_params * dyslider
controltop = panelbot + controlh
bmargin  = 0.15

fig, axs = plt.subplots(3, 1, figsize=(3.5, 3.5 * (3 + controlh) / 1.618), dpi=120)
axs = axs.flat
fig.subplots_adjust(left=0.2, top=0.95, bottom=controltop + bmargin, wspace=0.0, hspace=0.35)

# first axis: the size distributions

axs[0].set_xlabel('particle size [cm]')
axs[0].set_ylabel('$\sigma$ [g cm$^{-2}$]')
axs[0].set_ylim(1e-4, 1e-1)

# second axis: size-averaged opacities as fct of wavelength

axs[1].set_xlabel('wavelength [cm]')
axs[1].set_ylabel('$\kappa_\mathrm{abs}$ [cm$^2$/g]')

# third axis: mean opacities as function of temperature

axs[2].set_xlabel('$T$ [K]')
axs[2].set_ylabel('$\kappa_\mathrm{R/P}$ [cm$^2$/g]')
axs[2].set_ylim(1e-1, 1e4)

# NOTE: some lines are just place holders and the actual data will be plotted
# in the update function

# plot the filled areas

axs[2].fill_between(T_array, kP2min, kP2max, facecolor='C0', alpha=0.5)
axs[2].fill_between(T_array, kR2min, kR2max, facecolor='C1', alpha=0.5)

# plot the size distributions

line11, = axs[0].loglog(a, power_law, label='power-law')
line12, = axs[0].loglog(a, power_law, label='fit 1')
line13, = axs[0].loglog(a, power_law, label='fit 2')
line14, = axs[0].loglog([], [], 'k--')
leg1 = axs[0].legend(fontsize='xx-small')

# plot the size-averaged opacities

line21, = axs[1].loglog(lam, k_abs_p_w, label='power-law')
line22, = axs[1].loglog(lam, k_abs_p_w, label='fit 1')
line23, = axs[1].loglog(lam, k_abs_p_w, label='fit 2')
leg2 = axs[1].legend(fontsize='xx-small')

# fixed T_evap power-law mean opacities

line31, = axs[2].loglog(par_grid[0], k_P_pf, label='Planck')
line32, = axs[2].loglog(par_grid[0], k_R_pf, label='Rosseland')

dummy2, = axs[2].loglog([], [], '-', c='k', label='power-law')
dummy1, = axs[2].loglog([], [], '+', c='k', label='fit')
line33, = axs[2].loglog(par_grid[0][0], k_P_pf[0], '+', c=line31.get_color())
line34, = axs[2].loglog(par_grid[0][0], k_P_pf[0], '+', c=line31.get_color())
line35, = axs[2].loglog(par_grid[0][0], k_R_pf[0], '+', c=line32.get_color())
line36, = axs[2].loglog(par_grid[0][0], k_R_pf[0], '+', c=line32.get_color())
leg3 = axs[2].legend(fontsize='xx-small')

# ----------------
# MAKE THE SLIDERS
# ----------------

sliders = []
for ip in range(n_params):
    ax = fig.add_axes(
        [xslider, controltop - ip * dyslider, xslider + wslider, hslider],
        facecolor='lightgoldenrodyellow')

    slider = w.Slider(ax, par_names[ip], 0, n_grid[ip] - 1, valinit=par_grid[ip].searchsorted(par_ini[ip]), valfmt='%i')
    sliders.append(slider)

# add a button if coala is installed

if coala_installed:
    ax = ax = fig.add_axes(
        [xslider, controltop - n_params * dyslider, xslider + 0.2, hslider],
        facecolor='lightgoldenrodyellow')
    sim_button = w.Button(ax, 'simulate')

# -----------------
# CALLBACK FUNCTION
# -----------------


def update(val):

    indices = [int(s.val) for s in sliders]
    values = [grid[i] for grid, i in zip(par_grid, indices)]

    for slider, name, value in zip(sliders, par_names, values):
        slider.label.set_text(name + ' $ = {0:8.3g}$'.format(value))

    # update data

    res = get_data(values)

    # update plot

    line12.set_ydata(res['f1'])
    line13.set_ydata(res['f2'])

    line22.set_ydata(res['k_abs_f1'])
    line23.set_ydata(res['k_abs_f2'])

    line33.set_xdata(values[0])
    line34.set_xdata(values[0])
    line35.set_xdata(values[0])
    line36.set_xdata(values[0])

    line33.set_ydata(res['k_P_f1'])
    line34.set_ydata(res['k_P_f2'])
    line35.set_ydata(res['k_R_f1'])
    line36.set_ydata(res['k_R_f2'])


def update_sim(val):

    if val is False:
        return

    # get slider values

    indices = [int(s.val) for s in sliders]
    values = [grid[i] for grid, i in zip(par_grid, indices)]
    _T, _v_frag, _alpha, _sig_g = values

    # update data

    sim = coala.get_size_distri(_alpha, _T, _sig_g, d2g, r, _v_frag, M_star=M_star, rho_s=rho_s, amin=a[0], amax=a[-1])
    _sigd = sim['sig_d'][-1]
    _sigd = _sigd / _sigd.sum()
    line14.set_xdata(sim['a'])
    line14.set_ydata(_sigd)
    plt.draw()


# ------------------------------------------
# LINK CALLBACK AND AVOID GARBAGE COLLECTION
# ------------------------------------------

for s in sliders:
    s.on_changed(update)

fig._callback = [update]
fig._widgets = [sliders]

if coala_installed:
    sim_button.on_clicked(update_sim)
    fig._widgets.append(sim_button)

update(None)
update_sim(False)

plt.show()
