import numpy as np
import os

if not os.path.isdir('lnk'):
    os.mkdir('lnk')

lp = np.logspace(-5, 1, 200)

l_wn, n_wn = np.loadtxt('P94-waterice-n.csv').T
l_wk, k_wk = np.loadtxt('P94-waterice-k.csv').T
l_w = np.logspace(np.log10(max(l_wn[0], l_wk[0]) * 1e-4), np.log10(min(l_wn[-1], l_wk[-1]) * 1e-4), 100)
wn = np.interp(np.log10(l_w * 1e4), np.log10(l_wn), n_wn)
wk = np.interp(np.log10(l_w * 1e4), np.log10(l_wk), k_wk)
np.savetxt(os.path.join('lnk', 'P94-waterice.lnk'), np.vstack((l_w, wn, wk)).T)

l_tn, n_tn = np.loadtxt('P94-troilite-n.csv').T
l_tk, k_tk = np.loadtxt('P94-troilite-k.csv').T
l_t = np.logspace(np.log10(max(l_tn[0], l_tk[0]) * 1e-4), np.log10(min(l_tn[-1], l_tk[-1]) * 1e-4), 100)
tn = np.interp(np.log10(l_t * 1e4), np.log10(l_tn), n_tn)
tk = np.interp(np.log10(l_t * 1e4), np.log10(l_tk), k_tk)
np.savetxt(os.path.join('lnk', 'P94-troilite.lnk'), np.vstack((l_t, tn, tk)).T)

l_on, n_on = np.loadtxt('P94-organics-n.csv').T
l_ok, k_ok = np.loadtxt('P94-organics-k.csv').T
l_o = np.logspace(np.log10(max(l_on[0], l_ok[0]) * 1e-4), np.log10(min(l_on[-1], l_ok[-1]) * 1e-4), 100)
on = np.interp(np.log10(l_o * 1e4), np.log10(l_on), n_on)
ok = np.interp(np.log10(l_o * 1e4), np.log10(l_ok), k_ok)
np.savetxt(os.path.join('lnk', 'P94-organics.lnk'), np.vstack((l_o, on, ok)).T)
