import numpy as np
import os

if not os.path.isdir('lnk'):
    os.mkdir('lnk')

lp = np.logspace(-5, 1, 200)


for species in ['waterice', 'troilite', 'organics', 'olivine', 'iron', 'orthopyroxene']:

    l_n, n = np.loadtxt('P94-{}-n.csv'.format(species)).T
    l_k, k = np.loadtxt('P94-{}-k.csv'.format(species)).T
    l2 = np.logspace(np.log10(max(l_n[0], l_k[0]) * 1e-4), np.log10(min(l_n[-1], l_k[-1]) * 1e-4), 100)
    #
    # n interpolated in log-lin, except troilite and iron
    #
    if species in ['troilite', 'iron']:
        n2 = 10.**np.interp(np.log10(l2 * 1e4), np.log10(l_n), np.log10(n))
    else:
        n2 = np.interp(np.log10(l2 * 1e4), np.log10(l_n), n)
    #
    # k interpolated in log-log
    #
    k2 = 10.**np.interp(np.log10(l2 * 1e4), np.log10(l_k), np.log10(k))
    np.savetxt(os.path.join('lnk', 'P94-{}.lnk'.format(species)), np.vstack((l2, n2, k2)).T)
