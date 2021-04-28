import wsp_tools as wt
from wsp_tools import np, plt

T = 300000
print(np.isclose(wt.energy(T), 1.2993635678823885e-13))
print(np.isclose(wt.p(T), 3.365624812638601e-22))
print(np.isclose(wt.dB(T), 1.9687489006848795e-12))
print(np.isclose(wt.k(T), 3191460985702.0464))
print(np.isclose(wt.omega(T), 1.2321243049929178e21))
print(np.isclose(wt.v_g(T), 232796486.28087974))
print(np.isclose(wt.v_p(T), 386069048.1609881))
print(np.isclose(wt.zR(3191460985702.0464, 1e-6), 1.5957304928510232))
print(np.isclose(wt.w(0,1e-6,3191460985702.0464), 1e-6))
print(np.isclose(wt.roc(1,1e-6,3191460985702.0464), 3.5463558058145694))

X = np.linspace(-1e-8, 1e-8, 64)
Y = np.linspace(-1e-8, 1e-8, 128)
Z = np.linspace(-1e-8, 1e-8, 3)
x2, y2 = np.meshgrid(X, Y)
x3, y3, z3 = np.meshgrid(X, Y, Z)

bessel2 = wt.bessel(x2, y2, z = 0)
bessel3 = wt.bessel(x3, y3, z3)
besselPacket = wt.besselPacket(kres=128)

plt.imshow(np.abs(bessel2)**2)
plt.title("2d Bessel")
plt.show()

plt.imshow(np.abs(bessel3[:,:,0])**2)
plt.title("3d Bessel (slice)")
plt.show()

plt.imshow(np.abs(besselPacket[:,:,0])**2)
plt.title("Bessel Packet (slice)")
plt.show()

LG2 = wt.LG(x2, y2, z=0)
LG3 = wt.LG(x3, y3, z3)

plt.imshow(np.abs(LG2)**2)
plt.title("2d LG")
plt.show()

plt.imshow(np.abs(LG3[:,:,0])**2)
plt.title("3d LG (slice)")
plt.show()
