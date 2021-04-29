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

from . import constants as _
import numpy as np
from deprecated import deprecated

__all__ = ['abphase2d']

@deprecated(version='1.2.0', reason='You should now use wsp_tools.ltem.ab_phase.')
def abphase2d(mx, my, mz, Lx=1e-6, Ly=1e-6, p=np.array([0,0,1]), t=60e-9):
	"""DEPRECATED. You should instead use `wsp_tools.ltem.ab_phase.`

	Calculates the Aharonov-Bohm phase acquired by an electron passing through a 2d
	magnetization.

	**Parameters**

	* **mx** : _ndarray_ <br />
	The x-component of magnetization. Should be a 2-dimensional array.

	* **my** : _ndarray_ <br />
	The y-component of magnetization. Should be a 2-dimensional array.

	* **mz** : _ndarray_ <br />
	The z-component of magnetization. Should be a 2-dimensional array.

	* **Lx** : _number, optional_ <br />
	The x-length of the area calculated. (i.e., Lx = xmax - xmin). <br />
	Default is `Lx = 1e-6`.

	* **Ly** : _number, optional_ <br />
	The y-length of the area calculated. (i.e., Ly = ymax - ymin). <br />
	Default is `Ly = 1e-6`.

	* **p** : _ndarray, optional_ <br />
	A unit vector in the direction of electron propagation. Should be a 1darray with length 3. <br />
	Default is `p = np.array([0,0,1])` (the electron is propagating in the z-direction).

	* **t** : _number, optional_ <br />
	The thickness of the 2d magnetization. <br />
	Default is `t = 60e-9`.

	**Returns**

	* **phi** : _ndarray_
	The phase acquired by an electron passing through the specified magnetization.
	Output is a 2d array of the same shape as mx, my, and mz.
	"""
	Mx = np.fft.fft2(mx)
	My = np.fft.fft2(my)
	Mz = np.fft.fft2(mz)
	M = np.array([Mx, My, Mz])
	SI = np.fft.fftfreq(M[0].shape[1],Lx/M[0].shape[1])
	SJ = np.fft.fftfreq(M[0].shape[0],Ly/M[0].shape[0])
	si, sj = np.meshgrid(SI, SJ)
	s = np.array([si,sj,0*si])
	s_mag = np.sqrt(np.einsum('ijk,ijk->jk',s,s))
	s_mag[s_mag == 0] = np.max(s_mag)/s_mag.shape[0]/1000000

	sig = np.nan_to_num(s/s_mag, nan=0,posinf=0,neginf=0)
	Gts = 1/(np.einsum('i,ijk->jk',p,sig)**2 + p[2]**2) \
			* np.sinc(t*np.einsum('i,ijk->jk',p,sig)/p[2])
	weights = np.zeros_like(M[0], dtype=complex)

	pxpxM = np.cross(p, np.cross(p, M, axisb=0))
	sigxz = np.cross(sig, np.array([0,0,1]), axisa=0)
	d = np.einsum('ijk,ijk->ij', sigxz, pxpxM)
	weights = 1j * np.nan_to_num(t/s_mag, nan=0,neginf=0,posinf=0) * Gts * d

	phi = 2*_.e/_.hbar/_.c * np.fft.ifft2(weights)
	return(phi.real)

@deprecated(version='1.2.0', reason='You should now use wsp_tools.ltem.img_from_phase')
def propagate(cphase, dx = 1, dy = 1, defocus=0, wavelength=1.97e-12):
	"""Calculates the Lorentz image given a specific phase and defocus.

	This takes into account the microscope and phase transfer functions.

	**Parameters**

	* **x** : _ndarray_ <br />
	The x coordinates of the input complex phase. Dimension should be 2.

	* **y** : _ndarray_ <br />
	The y coordinates of the input complex phase. Dimension should be 2.

	* **cphase**: _complex ndarray_ <br />
	The complex phase to propagate. Dimension should be two. May be complex
	to allow for attenuation through the sample.

	* **defocus** : _number, optional_ <br />
	The defocus at which the phase was acquired. <br />
	Default is `defocus = 0`.

	* **wavelength** : _number, optional_ <br />
	The wavelength of the electron. Scales the phase transfer function and the reciprocal coordinates. <br />
	Default is `wavelength = 1.97e-12` (the relativistic wavelength of a 300kV electron).

	* **focal_length** : _number, optional_ <br />
	The focal length of the lens, which scales the reciprocal coordinates. <br />
	Default is `focal_length = 1`.

	**Returns**

	* **psi_out** : _complex ndarray_ <br />
	The transverse complex amplitude in the image plane. Output has the same
	shape as x, y, and cphase.
	"""
	U = np.fft.fftfreq(cphase.shape[1], dx)
	V = np.fft.fftfreq(cphase.shape[0], dy)
	qx, qy = np.meshgrid(U, V)
	psi_0 = np.exp(1j * cphase)
	psi_q = np.fft.fft2(psi_0)
	psi_out = np.fft.ifft2(psi_q * T(qx, qy, defocus, wavelength))
	return(psi_out)

def T(qx, qy, defocus, wavelength):
	"""Utility function for propagate(). Microscope transfer function.
	"""
	out = aperture(qx, qy) * np.exp(-1j * chi(qx, qy, defocus, wavelength))
	return(out)

def chi(qx, qy, defocus, wavelength):
	"""Utility function for propagate(). Phase transfer function.
	"""
	return(_.pi * wavelength * defocus * (qx**2 + qy**2))

def aperture(qx, qy):
	"""Utility function for propagate(). Circular aperture.
	"""
	out = np.zeros_like(qx)
	out[np.sqrt(qx**2 + qy**2) < np.max(np.sqrt(qx**2 + qy**2))] = 1
	return(out)

@deprecated(version='1.2.0', reason='You should now use wsp_tools.ltem.jchessmodel')
def jchessmodel(x, y, z=0, **kwargs):
	"""Calculates the magnetization of a hopfion based on Jordan Chess' model.

	**Parameters**

	* **x** : _number, ndarray_ <br />
	The x-coordinates over which to calculate magnetization.

	* **y** : _number, ndarray_ <br />
	The y-coordinates over which to calculate magnetization.

	* **z** : _number, ndarray, optional_ <br />
	The z-coordinates over which to calculate magnetization. Note, if z is an
	ndarray, then x, y, and z should have the same shape rather than relying
	on array broadcasting. <br />
	Default is `z = 0`.

	* **aa**, **ba**, **ca** : _number, optional_ <br />
	In this model, the thickness of the domain wall is set by a
	Gaussian function, defined as `aa * exp(-ba * z**2) + ca`. <br />
	Defaults are `aa = 5`, `ba = 5`, `ca = 0`.

	* **ak**, **bk**, **ck** : _number, optional_ <br />
	In this model, the thickness of the core is set by a Gaussian function,
	defined as `ak * exp(-bk * z**2) + ck`. <br />
	Defaults are `ak = 5e7`, `bk = -50`, `ck = 0`.

	* **bg**, **cg** : _number, optional_ <br />
	In this model, the helicity varies as a function of z, given
	as `pi / 2 * tanh( bg * z ) + cg`. <br />
	Defaults are `bg = 5e7`, `cg = pi/2`.

	* **n** : _number, optional_ <br />
	The skyrmion number. <br />
	Default is `n = 1`.


	**Returns**

	* **mx** : _ndarray_ <br />
	The x-component of magnetization. Shape will be the same as x and y.

	* **my** : _ndarray_ <br />
	The y-component of magnetization. Shape will be the same as x and y.

	* **mz** : _ndarray_ <br />
	The z-component of magnetization. Shape will be the same as x and y.
	"""
	p = {   'aa':5, 'ba':5, 'ca':0,
			'ak':5e7, 'bk':-5e1, 'ck':0,
			'bg':5e7, 'cg':_.pi/2, 'n': 1}
	for key in kwargs.keys():
		if not key in p.keys(): return("Error: {:} is not a kwarg.".format(key))
	p.update(kwargs)

	r, phi = np.sqrt(x**2+y**2), np.arctan2(y,x)

	alpha_z = p['aa'] * np.exp(-p['ba'] * z**2) + p['ca']
	k_z = p['ak'] * np.exp(-p['bk'] * z**2) + p['ck']
	gamma_z = _.pi / 2 * np.tanh(p['bg'] * z) + p['cg']
	Theta_rz = 2 * np.arctan2((k_z * r)**alpha_z, 1)

	mx = np.cos(p['n']*phi%(2*_.pi)-gamma_z) * np.sin(Theta_rz)
	my = np.sin(p['n']*phi%(2*_.pi)-gamma_z) * np.sin(Theta_rz)
	mz = np.cos(Theta_rz)
	return(mx, my, mz)
