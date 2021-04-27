import os, sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import wsp_tools as wt
from wsp_tools import np, plt

T = 300000
print(wt.energy(T))
print(wt.p(T))
print(wt.dB(T))
print(wt.k(T))
print(wt.omega(T))
print(wt.v_g(T))
print(wt.v_p(T))
print(wt.zR(wt.k(T), 1e-6))
print(wt.w(0,1e-6,wt.k(T)))
print(wt.roc(1,1e-6,wt.k(T)))

X = np.linspace(-1e-8, 1e-8, 128)
Z = np.linspace(-1e-8, 1e-8, 11)
x2, y2 = np.meshgrid(X, X)
x3, y3, z3 = np.meshgrid(X, X, Z)
bessel2 = wt.bessel(x2, y2, z = 0)
bessel3 = wt.bessel(x3, y3, z3)

plt.imshow(np.abs(bessel2)**2)
plt.title("2d Bessel")
plt.show()

plt.imshow(np.abs(bessel3[:,:,0])**2)
plt.title("3d Bessel")
plt.show()

besselPacket = wt.besselPacket()
plt.imshow(np.abs(besselPacket[:,:,0])**2)
plt.title("Bessel Packet")
plt.show()

LG2 = wt.LG(x2, y2, z=0)
LG3 = wt.LG(x3, y3, z3)

plt.imshow(np.abs(LG2)**2)
plt.title("2d LG")
plt.show()

plt.imshow(np.abs(LG3[:,:,0])**2)
plt.title("3d LG")
plt.show()
