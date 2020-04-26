# %%
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import gamma
from wsp_tools.constants import *
import wsp_tools.beam as beam
import gauge_invar_glaser as gig
plt.rcParams.update({'font.size':18})

# %%
######
fpath = '../images/'

######## 1. B0 = 0 - free solutions
w0 = 2e-6; k = beam.k(3e5); l=1

B0 = 0
a = .005
z0, zf = -14*a, 14*a
r = gig.solve_odes_LG(z0,zf,a=a, B0=B0, w0=w0, k=k, l=l)
fig, ((ax0, ax1), (ax2, ax3)) = plt.subplots(2,2,sharex='col',figsize=(15,10))
plt.suptitle("B(z) = 0")
ax0.plot(r.t, gig.Tp0f(r.t,w0,k,l),'--')
ax1.plot(r.t, gig.r02f(r.t,w0,k,l),'--')
ax2.plot(r.t, gig.Lz0f(r.t,w0,k,l),'--')
ax3.plot(r.t, gig.Gp0f(r.t,w0,k,l),'--')
ax0.plot(r.t, r.y[0],'-.')
ax1.plot(r.t, r.y[1],'-.')
ax2.plot(r.t, r.y[2],'-.')
ax3.plot(r.t, r.y[3],'-.')
ax0.set_title(r'$\left<T_{\perp}\right>$')
ax1.set_title(r'$\left<\rho^2\right>$')
ax2.set_title(r'$\left<L_z\right>$')
ax3.set_title(r'$\left<G_{\perp}\right>$')
plt.legend(['analytic','sim'])
plt.savefig(fpath+"B_0_analytic_and_sim.png")
plt.show()
