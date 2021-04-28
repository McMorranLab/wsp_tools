import wsp_tools as wt
from wsp_tools import np, plt

T = 300000

class TestBeamParams:
	def test_E(self):
		assert(np.isclose(wt.energy(T), 1.2993635678823885e-13))
	def test_p(self):
		assert(np.isclose(wt.p(T), 3.365624812638601e-22))
	def test_dB(self):
		assert(np.isclose(wt.dB(T), 1.9687489006848795e-12))
	def test_k(self):
		assert(np.isclose(wt.k(T), 3191460985702.0464))
	def test_omega(self):
		assert(np.isclose(wt.omega(T), 1.2321243049929178e21))
	def test_v_g(self):
		assert(np.isclose(wt.v_g(T), 232796486.28087974))
	def test_v_p(self):
		assert(np.isclose(wt.v_p(T), 386069048.1609881))
	def test_zR(self):
		assert(np.isclose(wt.zR(3191460985702.0464, 1e-6), 1.5957304928510232))
	def test_w(self):
		assert(np.isclose(wt.w(0,1e-6,3191460985702.0464), 1e-6))
	def test_roc(self):
		assert(np.isclose(wt.roc(1,1e-6,3191460985702.0464), 3.5463558058145694))

X = np.linspace(-1e-8, 1e-8, 64)
Y = np.linspace(-1e-8, 1e-8, 128)
Z = np.linspace(-1e-8, 1e-8, 3)
x2, y2 = np.meshgrid(X, Y)
x3, y3, z3 = np.meshgrid(X, Y, Z)

class TestBesselShape:
	def test_bessel2(self):
		bessel2 = wt.bessel(x2, y2, z = 0)
		assert(bessel2.shape == (128, 64))
	def test_bessel3(self):
		bessel3 = wt.bessel(x3, y3, z3)
		assert(bessel3.shape == (128, 64, 3))

class TestBesselPacketShape:
	def test_bessel_packet(self):
		besselPacket = wt.besselPacket(kres=4)
		assert(besselPacket.shape == (4, 4, 4))

class TestLGShape:
	def test_LG2(self):
		LG2 = wt.LG(x2, y2, z=0)
		assert(LG2.shape == (128, 64))
	def test_LG3(self):
		LG3 = wt.LG(x3, y3, z3)
		assert(LG3.shape == (128, 64, 3))
