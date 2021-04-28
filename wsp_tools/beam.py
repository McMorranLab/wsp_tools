# wsp-tools is TEM data analysis and simulation tools developed by WSP as a grad student in the McMorran Lab.
# Copyright (C) 2021  William S. Parker
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

"""Module to generate spatial modes and calculate beam parameters.
"""

import numpy as np
from . import constants as _
from scipy.special import eval_genlaguerre, factorial

__all__ = ['energy','p','dB','k','omega','v_g','v_p','bessel','besselPacket','zR','w','roc','LG']

def energy(T_eV):
	"""Calculates the total electron energy from its kinetic energy.

	**Parameters**

	* **T_eV** : _number_ <br />
	the kinetic energy in eV.

	**Returns**

	* **energy(T_eV)** : _number_ <br />
	the total energy in eV.
	"""
	return(T_eV * _.e + _.m_e * _.c**2)

def p(T_eV):
	"""Calculates the momentum from the electron's kinetic energy.

	**Parameters**

	* **T_eV** : _number_ <br />
	the kinetic energy in eV.

	**Returns**

	* **p(T_eV)** : _number_ <br />
	the momentum in N s.
	"""
	return(1/_.c * np.sqrt(energy(T_eV)**2 - (_.m_e * _.c**2)**2))

def dB(T_eV):
	"""Calculates the de Broglie wavelength from the electron's kinetic energy.

	**Parameters**

	* **T_eV** : _number_ <br />
	the kinetic energy in eV.

	**Returns**

	* **dB(T_eV)** : _number_ <br />
	the de Broglie wavelength in m.
	"""
	return(_.h / p(T_eV))

def k(T_eV):
	"""Calculates the wavenumber from the electron's kinetic energy.

	**Parameters**

	* **T_eV** : _number_ <br />
	the kinetic energy in eV.

	**Returns**

	* **k(T_eV)** : _number_ <br />
	the wavenumber in m<sup>-1</sup>.
	"""
	return(p(T_eV) / _.hbar)

def omega(T_eV):
	"""Calculates the angular frequency from the electron's kinetic energy.

	**Parameters**

	* **T_eV** : _number_ <br />
	the kinetic energy in eV.

	**Returns**

	* **omega(T_eV)** : _number_ <br />
	the angular frequency in rad s<sup>-1</sup>.
	"""
	return(energy(T_eV) / _.hbar)

def v_g(T_eV):
	"""Calculates the group velocity from the electron's kinetic energy.

	**Parameters**

	* **T_eV** : _number_ <br />
	the kinetic energy in eV.

	**Returns**

	* **v_g(T_eV)** : _number <br />
	the group velocity in m s<sup>-1</sup>.
	"""
	return(_.c * k(T_eV) / np.sqrt(_.m_e**2 * _.c**2 / _.hbar**2 + k(T_eV)**2))

def v_p(T_eV):
	"""Calculates the phase velocity from the electron's kinetic energy.

	**Parameters**

	* **T_eV** : _number_ <br />
	the kinetic energy in eV.

	**Returns**

	* **v_p(T_eV)** : _number_ <br />
	the phase velocity in m s<sup>-1</sup>.
	"""
	return(omega(T_eV)/k(T_eV))

def bessel(x, y, z = 0, l = 0, kz0 = k(3e5),
			kperp = 0.0005*k(3e5), dkperp = 0, N = 30):
	"""Creates a Bessel beam by adding plane waves.

	The spectrum is a circle in k-space <br />`(kz0, dkperp + kperp * cos(theta), kperp * sin(theta))`.

	**Parameters**

	* **x** : _number, ndarray_ <br />
	the x-coordinates over which to calculate the beam.

	* **y** : _number, ndarray_ <br />
	the y-coordinates over which to calculate the beam.

	* **z** : _number, ndarray, optional_ <br />
	the z-coordinates over which to calculate the beam. <br />
	Default is `z = 0`.

	* **l** : _number, optional_ <br />
	the winding number, or chirality, of the beam. <br />
	Default is `l = 0`.

	* **kz0** : _number, optional_ <br />
	the z-coordinate of the spectrum. That is, the z component of the beam's wavevector. <br />
	Default is `kz0 = k(3e5)`.

	* **kperp** : _number, optional_ <br />
	the perpendicular component of the spectrum. That is, the perpendicular component of the beam's wavevector. This is the radius of the beam's spectrum. <br />
	Default is `kperp = 5e-4 * k(3e5)`.

	* **dkperp** : _number, optional_ <br />
	the perpendicular offset of the spectrum. This is applied in the `k_x` direction. <br />
	Default is `dkperp = 0`.

	* **N** : _int, optional_ <br />
	The number of plane waves to use. <br />
	Default is `N = 30`.

	**Returns**

	* **mode** : _ndarray_ <br />
	The complex amplitudes of the Bessel beam at the specified x, y, and z positions.
	Shape is broadcast between x, y, and z inputs.
	"""
	mode = x*y*z*0j
	w_0 = np.sqrt((kz0**2 + kperp**2)*_.c**2 + _.m_e**2*_.c**4/_.hbar**2)
	for n in range(N):
		phi_n = 2*np.pi*n/N
		kxn, kyn, kzn = kperp*np.cos(phi_n) + dkperp, kperp*np.sin(phi_n), kz0
		w_n = np.sqrt(w_0**2 + dkperp**2 * _.c**2 + 2*kxn*dkperp*_.c**2)
		mode += np.exp(1j*(kxn*x + kyn*y + kzn*z + l*phi_n))
	mode *= 1/np.sqrt(N)
	return(mode)

def besselPacket(t = 0, l = 0,
		kres = 2**7, kmin = -3*k(3e5), kmax = 3*k(3e5),
		kz0 = k(3e5), kperp = .5*k(3e5), dkperp = 0,
		sig = 0.05*k(3e5)):
	"""Creates a bessel beam by Fourier transforming a Gaussian spectrum.

	The spectrum is a gaussian centered on a circle in
	k-space <br /> `(k_z, dkperp + kperp * cos(theta), kperp * sin(theta))`.

	**Parameters**

	* **t** : _number, optional_ <br />
	time in seconds. <br />
	Default is `t = 0`.

	* **l** : _number, optional_ <br />
	the winding number, or chirality, of the beam. <br />
	Default is `l = 0`.

	* **kres** : _int, optional_ <br />
	the resolution of the k-space; also, the resolution of the output beam. Note that the function will be faster if the resolution is a power of 2. <br />
	Default is `kres = 128`.

	* **kmin** : _number, optional_ <br />
	the minimum value in k-space. This is applied to the x, y, and z components of the wavenumber. <br />
	Default is `kmin = -3 * k(3e5)`.

	* **kmax** : _number, optional_ <br />
	the maximum value in k-space. This is applied to the x, y, and z components of the wavenumber. <br />
	Default is `kmax = 3 * k(3e5)`.

	* **kz0** : _number, optional_ <br />
	the z-coordinate of the spectrum. That is, the z component of the beam's wavevector. <br />
	Default is `kz0 = k(3e5)`.

	* **kperp** : _number, optional_ <br />
	the perpendicular component of the spectrum. That is, the perpendicular component of the beam's wavevector. This is the radius of the beam's spectrum. <br />
	Default is `kperp = 5e-4 * k(3e5)`.

	* **dkperp** : _number, optional_ <br />
	the perpendicular offset of the spectrum. This is applied in the `k_x` direction. <br />
	Default is `dkperp = 0`.

	* **sig** : _number, optional_ <br />
	the standard deviation of the Gaussian envelope around the circle in k-space. <br />
	Default is `sig = 0.05 * k(3e5)`.

	**Returns**

	* **mode** : _ndarray_ <br />
	three-dimensional array containing the complex amplitudes of the Bessel packet. Shape is `kres`x`kres`x`kres`.
	"""
	# K = np.linspace(kmin, kmax, kres)
	kx, ky, kz = np.ogrid[kmin:kmax:1j*kres,kmin:kmax:1j*kres,kmin:kmax:1j*kres]
	w = _.c * np.sqrt((kx**2 + ky**2 + kz**2) + _.c**2 * _.m_e**2 /_.hbar**2)
	# cylindrical
	kr = np.sqrt(kx**2 + ky**2)
	mode = kz*0j
	phi = np.arctan2(ky,kx)
	spec = np.exp( -( (kr-kperp)**2 + (kz-kz0)**2 ) / (2*sig**2)) \
			* np.exp(1j*phi*l-1j*w*t)
	### Include first ifftshift because spec is [-k0, k0], not [0,k0]
	mode = np.fft.fftshift(np.fft.ifftn(np.fft.ifftshift(spec)))
	return(mode)

def zR(k, w0):
	"""Rayleigh range as a function of wavenumber and beam waist.

	**Parameters**

	* **k** : _number_ <br />
	Wavenumber.

	* **w0** : _number_ <br />
	Beam waist.

	**Returns**

	* **zR(k, w0)** : _number_ <br />
	the Rayleigh range.
	"""
	lam = 2*_.pi/k
	return(_.pi * w0**2 / lam)

def w(z,w0,k):
	"""Spot size as a function of z, beam waist, and wavenumber.

	**Parameters**

	* **z** : _number_ <br />
	z-position.

	* **w0** : _number_ <br />
	beam waist.

	* **k** : _number_ <br />
	wavenumber.

	**Returns**

	* **w(z, w0, k)** : _number_ <br />
	the spot size at z.
	"""
	return(w0 * np.sqrt(1 + (z/zR(k,w0))**2))

np.seterr(divide='ignore')
def roc(z, w0, k):
	"""Radius of curvature as a function of z, beam waist, and wavenumber.

	**Parameters**

	* **z** : _number_ <br />
	z-position.

	* **w0** : _number_ <br />
	beam waist.

	* **k** : _number_ <br />
	wavenumber.

	**Returns**

	* **roc(z, w0, k)** : _number_ <br />
	Radius of curvature.
	"""
	return(z + np.divide(zR(k,w0)**2,z))

def LG(x, y, z = 0, l = 0, p = 0, w_0 = 2e-6, lam = 1.97e-12):
	"""Generates a Laguerre-Gauss spatial mode.

	Note: broadcasting is not supported - if x and y are both 1d arrays, the
	result will be a 1darray and will NOT make sense. Likewise if x, y are 2d
	arrays, but z is a 1d array, the result will be an error.

	**Parameters**

	* **x** : _number, ndarray_ <br />
	The x-coordinates over which to calculate the beam.

	* **y** : _number, ndarray_ <br />
	The y-coordinates over which to calculate the beam.

	* **z** : _number, ndarray, optional_ <br />
	The z-coordinates over which to calculate the beam. <br />
	Default is `z = 0`.

	* **l** : _number, optional_ <br />
	The winding number, or chirality, of the beam. <br />
	Default is `l = 0`.

	* **p** : _number, optional_ <br />
	The radial index of the beam. <br />
	Default is `p = 0`.

	* **w_0** : _number, optional_ <br />
	The beam waist. <br />
	Default is `w_0 = 2e-6`.

	* **lam** : _number, optional_ <br />
	The beam's wavelength. <br />
	Default is `lam = 1.97e-12` (the relativistic wavelength of a 300keV electron).

	**Returns**

	* **mode** : _ndarray_ <br />
	The complex amplitudes of a Laguerre Gaussian beam at the specified x, y, and z positions.
	"""
	r, theta = np.sqrt(x**2 + y**2), np.arctan2(y,x)
	Clp = np.sqrt(2*factorial(p)/np.pi/factorial(p+np.abs(l)))
	z_r = np.pi*w_0**2/lam
	w_z = w_0*np.sqrt(1 + z**2/z_r**2)
	gouy = (2*p + np.abs(l) + 1) * np.nan_to_num(np.arctan2(z,z_r))
	r = np.sqrt(x**2+y**2); theta = np.arctan2(y,x)
	mode = Clp \
			* w_0/w_z \
			* (r * np.sqrt(2)/w_z)**np.abs(l) \
			* np.exp(-r**2/w_z**2) \
			* eval_genlaguerre(p,np.abs(l),(2*r**2)/w_z**2) \
			* np.exp(-1j * 2 * np.pi / lam * r**2 * z/ 2 /(z**2 + z_r**2)) \
			* np.exp(-1j * l * theta) \
			* np.exp(1j * (1 + np.abs(l) * 2*p))
	mode /= np.max(np.abs(mode))
	return(mode)
