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
##### Glaser potential - 2x2 plotting
fpath = "../images/image_seq_hor/"

#### Beam parameters
w0 = 2e-6; k = beam.k(3e5); l=1

#### B-field parameters
a = .005
z0, zf = -14*a, 20*a
Bs = np.linspace(.01,3,300,endpoint=True)

####### Plotting command
n = 0
for B0 in Bs:
    n += 1
    r = gig.solve_odes_LG(z0,zf,a=a, B0=B0, w0=w0, k=k, l=l)
    fig, ((ax0,ax1),(ax2,ax3)) = plt.subplots(2,2,sharex='col',figsize=(15,10))
    sel = np.abs(r.y[1])/r.y[1][0] < 1.5
    plt.suptitle(r"$B_0$ = {:.4f} T".format(B0),fontsize="32")
    # ax0.plot(r.t/a, np.zeros_like(r.t),'black')
    # ax1.plot(r.t/a, np.zeros_like(r.t),'black')
    # ax2.plot(r.t/a, np.zeros_like(r.t),'black')
    # ax3.plot(r.t/a, np.zeros_like(r.t),'black')
    ax0.plot(r.t/a, r.y[0], 'tab:blue')
    ax1.plot(r.t[sel]/a, r.y[1][sel], 'tab:green')
    ax2.plot(r.t/a, r.y[2], 'tab:orange')
    ax3.plot(r.t/a, r.y[3], 'tab:red')
    ax0.set_ylim([-.2e-23,3.2e-23])
    ax1.set_ylim([-7e-12,7e-12])
    ax2.set_ylim([-1.6e-31,1.5e-31])
    ax3.set_ylim([-1.6e-31,1.5e-31])
    ax0.plot(r.t/a, np.zeros_like(r.t), 'black', linewidth=.3)
    ax0.plot(np.zeros_like(r.t),
            np.linspace(ax0.get_ylim()[0],ax0.get_ylim()[1],len(r.t)),
            'black', linewidth=.3)
    ax1.plot(r.t/a, np.zeros_like(r.t), 'black', linewidth=.3)
    ax1.plot(np.zeros_like(r.t),
            np.linspace(ax1.get_ylim()[0],ax1.get_ylim()[1],len(r.t)),
            'black', linewidth=.3)
    ax2.plot(r.t/a, np.zeros_like(r.t), 'black', linewidth=.3)
    ax2.plot(np.zeros_like(r.t),
            np.linspace(ax2.get_ylim()[0],ax2.get_ylim()[1],len(r.t)),
            'black', linewidth=.3)
    ax3.plot(r.t/a, np.zeros_like(r.t), 'black', linewidth=.3)
    ax3.plot(np.zeros_like(r.t),
            np.linspace(ax3.get_ylim()[0],ax3.get_ylim()[1],len(r.t)),
            'black', linewidth=.3)
    ax0.set_title(r'$\left<T_{\perp}\right>$')
    ax1.set_title(r'$\left<\rho^2\right>$')
    ax2.set_title(r'$\left<L_z\right>$')
    ax3.set_title(r'$\left<G_{\perp}\right>$')
    ax2.set_xlabel('z/a'); ax3.set_xlabel('z/a')
    plt.show()
    # plt.savefig(fpath+"img_{:04d}.png".format(n))
    plt.close(fig)
    # print(n)
