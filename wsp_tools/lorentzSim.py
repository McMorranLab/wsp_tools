from . import constants as _
from . import np

__all__ = ['ab_phase','abphase2d','propagate','T','chi','aperture','jchessmodel','A_from_mag','B_from_mag']

def G(p, sig, z_hat, ts_mag):
	sum1 = np.einsum('i,i...->...',p,sig)
	sum2 = np.einsum('i,i...->...', p, z_hat)
	out = 1 / (sum1**2 + sum2**2)[np.newaxis,...]
	out *= np.sinc(ts_mag * sum1 / sum2)
	return(out)

def ab_phase(mx, my, mz, dx=1, dy=1, thickness=60e-9, p = np.array([0,0,1])):
	"""Calculates the Aharanov-Bohm phase shift of an electron passing through a material with specified magnification.

	This is an implementation of the derivation presented by Mansuripur et al, 1991, _Computation of electron diffraction patterns in Lorentz TEM_.

	The original derivation is for a single slice of material, with magnetization constant through its thickness, but varying in x and y. This implementation can perform that calculation on multiple slices at once (i.e., mx, my, and mz can be 3d arrays, rather than 2d). In this case, the output will have the same shape as mx, my, and mz, such that the total phase shift will be the sum of the output along its last axis.

	If the electron is incident on the material at an oblique angle, **this should be reflected in the input magnetization**. That is, this implementation does **not** take into account that the electron will pass through each slice of material with an offset in x and/or y. The unit vector p describing the electron's motion is only used for each slice.

	**Parameters**

	* **mx** : _ndarray_ <br />
	A 2d array describing the x-component of magnetization in the film.

	* **my** : _ndarray_ <br />
	A 2d array describing the y-component of magnetization in the film.

	* **mz** : _ndarray_ <br />
	A 2d array describing the z-component of magnetization in the film.

	* **dx** : _scalar_ <br />
	The spacing between x-values.

	* **dy** : _scalar_ <br />
	The spacing between y-values.

	* **thickness** : _scalar_ <br />
	The thickness of each magnetization slice.

	* **p** : _ndarray_ <br />
	A 1d array of length 3, constituting a unit vector in the direction of the electron's motion.
	"""
	mx = np.atleast_3d(mx)
	my = np.atleast_3d(my)
	mz = np.atleast_3d(mz)
	Mx = np.fft.fft2(mx, axes=(0,1))
	My = np.fft.fft2(my, axes=(0,1))
	Mz = np.fft.fft2(mz, axes=(0,1))
	M = np.array([Mx, My, Mz])

	Sx = np.fft.fftfreq(mx.shape[1], dx)
	Sy = np.fft.fftfreq(mx.shape[0], dy)
	sx, sy = np.meshgrid(Sx, Sy)

	s = np.array([sx, sy, 0*sy])[...,np.newaxis]
	s_mag = np.sqrt(np.einsum('i...,i...->...',s,s))[np.newaxis,...]
	sig = s/s_mag

	z_hat = np.array([np.zeros_like(mx), np.zeros_like(mx), np.ones_like(mx)])

	Gp = G(p, sig, z_hat, thickness * s_mag)
	sig_x_z = np.cross(sig, z_hat, axisa=0, axisb=0, axisc=0)
	p_x_p_M = np.cross(p, np.cross(p, M, axisa=0, axisb=0, axisc=0), axis=0, axisb=0, axisc=0)
	weights = 1j * thickness / s_mag * Gp * np.einsum('i...,i...->...',sig_x_z, p_x_p_M)
	weights[:,0,0,:] = 0
	phase = 2 * _.e / _.hbar / _.c * np.fft.ifft2(weights, axes=(1,2))
	return(np.squeeze(phase.real))

def A_from_mag(
		mx, my, mz, z = 0, dx = 1, dy = 1,
		thickness = 60e-9
		):
	"""Calculates the magnetic vector potential for a given 2d magnetization, at a particular z.

	The calculation is an implementation of the math described in Mansuripur et. al, 1991 - _Computation of electron diffraction patterns in Lorentz TEM_.

	Note: the result is an array of shape (3, mx.shape[0], mx.shape[1], len(z)) - this can get very large when many z values are included, and it can be more efficient to run a for loop over z than to generate an array that takes up that much memory.

	**Parameters**

	* **mx** : _ndarray_ <br />
	A 2d array describing the x-component of magnetization in the film.

	* **my** : _ndarray_ <br />
	A 2d array describing the y-component of magnetization in the film.

	* **mz** : _ndarray_ <br />
	A 2d array describing the z-component of magnetization in the film.

	* **z** : _scalar, ndarray_ <br />
	A scalar or 1d array describing the z-value at which to calculate the magnetic vector potential.

	* **dx** : _scalar_ <br />
	The spacing between x-values.

	* **dy** : _scalar_ <br />
	The spacing between y-values.

	* **thickness** : _scalar_ <br />
	The thickness of the magnetization slice.

	**Returns**

	* **A**: _ndarray_ <br />
	The magnetic vector potential. If **z** is a scalar, its shape will be (3, mx.shape[0], mx.shape[1]).
	If **z** is a 1d array, its shape will be (3, mx.shape[0], mx.shape[1], len(z)).
	"""
	z = np.atleast_1d(z)
	selz_m = z < -thickness/2
	selz_z = np.abs(z) <= thickness/2
	selz_p = z > thickness/2
	z = z[np.newaxis,np.newaxis,np.newaxis,...]

	# All the np.newaxis stuff is so that everything has the following axes:
	# (3 vector components, y, x, z) and is therefore broadcastable

	z_hat = np.array([np.zeros_like(mx), np.zeros_like(mx), np.ones_like(mx)])[...,np.newaxis]

	Mx = np.fft.fft2(mx)
	My = np.fft.fft2(my)
	Mz = np.fft.fft2(mz)
	M = np.array([Mx, My, Mz])[...,np.newaxis]

	Sx = np.fft.fftfreq(mx.shape[1], dx)
	Sy = np.fft.fftfreq(mx.shape[0], dy)
	sx, sy = np.meshgrid(Sx, Sy)

	s = np.array([sx, sy, 0*sy])[...,np.newaxis]
	s_mag = np.sqrt(np.einsum('i...,i...->...',s,s))[np.newaxis,...]
	sig = s/s_mag
	sigp = sig + 1j * z_hat
	sigm = sig - 1j * z_hat

	A_mn = get_A_mn_components(
				mx.shape[0], mx.shape[1], z.shape[-1], selz_m, selz_z,
				selz_p, s_mag, z, sigm, sigp, sig, M, thickness, z_hat)
	A = np.fft.ifft2(A_mn, axes=(1,2))
	return(np.squeeze(A))

def B_from_mag(
		mx, my, mz, z = 0, dx = 1, dy = 1,
		thickness = 60e-9
		):
	"""Calculates the magnetic field for a given 2d magnetization, at a particular z.

	The calculation is an implementation of the math described in Mansuripur et. al, 1991 - _Computation of electron diffraction patterns in Lorentz TEM_.

	Note: the result is an array of shape (3, mx.shape[0], mx.shape[1], len(z)) - this can get very large when many z values are included, and it can be more efficient to run a for loop over z than to generate an array that takes up that much memory.

	**Parameters**

	* **mx** : _ndarray_ <br />
	A 2d array describing the x-component of magnetization in the film.

	* **my** : _ndarray_ <br />
	A 2d array describing the y-component of magnetization in the film.

	* **mz** : _ndarray_ <br />
	A 2d array describing the z-component of magnetization in the film.

	* **z** : _scalar, ndarray_ <br />
	A scalar or 1d array describing the z-value at which to calculate the magnetic vector potential.

	* **dx** : _scalar_ <br />
	The spacing between x-values.

	* **dy** : _scalar_ <br />
	The spacing between y-values.

	* **thickness** : _scalar_ <br />
	The thickness of the magnetization slice.

	**Returns**

	* **B**: _ndarray_ <br />
	The magnetic field. If **z** is a scalar, its shape will be (3, mx.shape[0], mx.shape[1]).
	If **z** is a 1d array, its shape will be (3, mx.shape[0], mx.shape[1], len(z)).
	"""
	### calculates the magnetic vector potential for a single slice of material
	### mx, my, mz are 2d arrays
	### z is a scalar or 1d array of the z positions at which to calculate A
	z = np.atleast_1d(z)
	selz_m = z < -thickness/2
	selz_z = np.abs(z) <= thickness/2
	selz_p = z > thickness/2
	z = z[np.newaxis,np.newaxis,np.newaxis,...]

	Mx = np.fft.fft2(mx)
	My = np.fft.fft2(my)
	Mz = np.fft.fft2(mz)
	M = np.array([Mx, My, Mz])[...,np.newaxis]

	Sx = np.fft.fftfreq(mx.shape[1], dx)
	Sy = np.fft.fftfreq(mx.shape[0], dy)
	sx, sy = np.meshgrid(Sx, Sy)

	s = np.array([sx, sy, 0*sy])[...,np.newaxis]
	s_mag = np.sqrt(np.einsum('i...,i...->...',s,s))[np.newaxis,...]
	sig = s/s_mag

	z_hat = np.array([np.zeros_like(mx), np.zeros_like(mx), np.ones_like(mx)])[...,np.newaxis]
	sigp = sig + 1j * z_hat
	sigm = sig - 1j * z_hat

	A_mn = get_A_mn_components(
			mx.shape[0], mx.shape[1], z.shape[-1], selz_m, selz_z,
			selz_p, s_mag, z, sigm, sigp, sig, M, thickness, z_hat)

	dxA_mn = 1j * 2 * _.pi * s[0][np.newaxis, ...] * A_mn
	dyA_mn = 1j * 2 * _.pi * s[1][np.newaxis, ...] * A_mn
	dzA_mn = np.zeros((3, mx.shape[0], mx.shape[1], z.shape[-1]), dtype=complex)
	dzA_mn[...,selz_m] = 2 * _.pi * s_mag * A_mn[...,selz_m]
	dzA_mn[...,selz_p] = -2 * _.pi * s_mag * A_mn[...,selz_p]
	dzA_mn[...,selz_z] = (2 * 1j / s_mag * np.cross((
						sig
						- 0.5 * 2 * _.pi * s_mag * np.exp(2 * _.pi * s_mag * (z[...,selz_z] - thickness / 2)) * sigm
						+ 0.5 * 2 * _.pi * s_mag * np.exp(-2 * _.pi * s_mag * (z[...,selz_z] + thickness / 2)) * sigp
						), M, axisa=0, axisb=0, axisc=0))
	dzA_mn[:,0,0,selz_m] = 0
	dzA_mn[:,0,0,selz_p] = 0
	dzA_mn[:,0,0,selz_z] = -4 * _.pi
	B_mn = np.zeros((3, mx.shape[0], mx.shape[1], z.shape[-1]), dtype=complex)
	B_mn[0] = dyA_mn[2] - dzA_mn[1]
	B_mn[1] = dzA_mn[0] - dxA_mn[2]
	B_mn[2] = dxA_mn[1] - dyA_mn[0]
	B = np.fft.ifft2(B_mn, axes=(1,2))
	return(np.squeeze(B))

def get_A_mn_components(
		xshape, yshape, zshape, selz_m, selz_z,
		selz_p, s_mag, z, sigm, sigp, sig, M, thickness, z_hat):
	### set everything to be broadcastable so we can write the A_mn equations
	### shape is: (3 vector components, y-axis, x-axis, z-axis) (thank meshgrid for switching x and y)

	A_mn = np.zeros((3, xshape, yshape, zshape), dtype=complex)
	A_mn[...,selz_m] = (2 * 1j / s_mag
						* np.exp(2 * _.pi * s_mag * z[...,selz_m])
						* np.sinh(_.pi * thickness * s_mag)
						* np.cross(sigm, M, axisa=0, axisb=0, axisc=0))
	A_mn[...,selz_z] = (2 * 1j / s_mag * np.cross((
						sig
						- 0.5 * np.exp(2 * _.pi * s_mag * (z[...,selz_z] - thickness / 2)) * sigm
						- 0.5 * np.exp(-2 * _.pi * s_mag * (z[...,selz_z] + thickness / 2)) * sigp
						), M, axisa=0, axisb=0, axisc=0))
	A_mn[...,selz_p] = (2 * 1j / s_mag
						* np.exp(-2 * _.pi * s_mag * z[...,selz_p])
						* np.sinh(_.pi * thickness * s_mag)
						* np.cross(sigp, M, axisa=0, axisb=0, axisc=0)
						)
	zero_comp = np.cross(z_hat[:,0,0,:], M[:,0,0,:], axisa=0, axisb=0, axisc=0)
	A_mn[:,0,0,selz_z] = -4 * _.pi * z[:,0,0,selz_z] * zero_comp
	A_mn[:,0,0,selz_m] = 2 * _.pi * thickness * zero_comp
	A_mn[:,0,0,selz_p] = -2 * _.pi * thickness * zero_comp
	return(A_mn)

def abphase2d(mx, my, mz, Lx=1e-6, Ly=1e-6, p=np.array([0,0,1]), t=60e-9):
	"""Calculates the Aharonov-Bohm phase acquired by an electron passing through a 2d
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

def propagate(x, y, cphase, defocus=0, wavelength=1.97e-12, focal_length=1):
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
	dx = x[1,1]-x[0,0]
	dy = y[1,1]-y[0,0]
	U = np.fft.fftfreq(cphase.shape[1],dx)
	V = np.fft.fftfreq(cphase.shape[0],dy)
	qx, qy = np.meshgrid(U, V)
	psi_0 = np.exp(1j*cphase)
	psi_q = np.fft.fft2(psi_0)
	psi_out = np.fft.ifft2(psi_q * T(qx, qy, defocus, wavelength))
	return(psi_out)

def T(qx, qy, defocus, wavelength):
	"""Utility function for propagate(). Microscope transfer function.
	"""
	out = aperture(qx, qy) * np.exp(-1j*chi(qx, qy, defocus, wavelength))
	return(out)

def chi(qx, qy, defocus, wavelength):
	"""Utility function for propagate(). Phase transfer function.
	"""
	return(_.pi*wavelength*defocus*(qx**2 + qy**2))

def aperture(qx, qy):
	"""Utility function for propagate(). Circular aperture.
	"""
	out = np.zeros_like(qx)
	out[np.sqrt(qx**2 + qy**2) < np.max(np.sqrt(qx**2 + qy**2))] = 1
	return(out)

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
