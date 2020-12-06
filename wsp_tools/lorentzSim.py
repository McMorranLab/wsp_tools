from . import constants as _
from . import np

def abphase2d(mx, my, mz, Lx=1e-6, Ly=1e-6, p=np.array([0,0,1]), t=60e-9):
	"""Calculates the Aharonov-Bohm phase acquired by an electron through a 2d
	magnetization.

	Inputs:

	The real space magnetization of the sample.

	mx: numpy.2darray() - x component of magnetization

	my: numpy.2darray() - y component of magnetization

	mz: numpy.2darray() - z component of magnetization

	Lx: scalar, the width of the sample.

	Ly: scalar, the length of the sample.

	p: numpy.1darray with length 3 - unit vector in the direction of electron motion.

	t: the thickness of the sample.

	Returns:

	phi = numpy.2darray() the phase acquired.
	"""
	Mx = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(mx)))
	My = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(my)))
	Mz = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(mz)))
	M = np.array([Mx, My, Mz])
	xres, yres = int(M.shape[1]/2), int(M.shape[2]/2)
	SI = np.arange(-xres,xres)
	SJ = np.arange(-yres,yres)
	si, sj = np.meshgrid(SI, SI)
	s = np.array([si/Lx,sj/Ly,0*si])
	s_mag = np.sqrt(np.einsum('ijk,ijk->jk',s,s))
	sig = np.nan_to_num(s/s_mag, nan=0,posinf=0,neginf=0)
	Gts = 1/(np.einsum('i,ijk->jk',p,sig)**2 + p[2]**2) \
			* np.sinc(t*np.einsum('i,ijk->jk',p,sig)/p[2])
	weights = np.zeros_like(M[0], dtype=complex)

	pxpxM = np.cross(p, np.cross(p, M, axisb=0))
	sigxz = np.cross(sig, np.array([0,0,1]), axisa=0)
	d = np.einsum('ijk,ijk->ij', sigxz, pxpxM)
	weights = 1j * np.nan_to_num(t/s_mag, nan=0,neginf=0,posinf=0) * Gts * d

	phi = 2*_.e/_.hbar/_.c * np.fft.fftshift(np.fft.ifft2(np.fft.fftshift(weights)))
	return(np.real(phi))

def propagate(x, y, cphase, defocus=0, wavelength=1.97e-12, focal_length=1):
	"""Calculates the Lorentz image given a specific phase.

	Inputs:

	x, y: numpy.2darray objects with the x and y coordinates

	cphase: the complex phase acquired (complex to allow for attenuation)

	defocus: scalar

	Returns:

	psi_out: numpy.2darray containing the complex wavefunction in the image plane.
	"""
	dx = x[1,1]-x[0,0]
	dy = y[1,1]-y[0,0]
	U = np.fft.fftshift(np.fft.fftfreq(x.shape[0],dx))
	V = np.fft.fftshift(np.fft.fftfreq(x.shape[1],dy))
	qx, qy = np.meshgrid(U, V)
	psi_0 = np.exp(1j*cphase)
	psi_q = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(psi_0)))
	psi_out = np.fft.fftshift(np.fft.ifft2(np.fft.fftshift(psi_q * T(qx, qy, defocus, wavelength))))
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

	Returns:

	mx, my, mz - numpy.ndarray objects in the same shape as x and y.
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
