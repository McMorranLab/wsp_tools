from numpy.fft import fftshift, fft2, ifft2, fftfreq
from numpy import exp, array, arange, meshgrid, sqrt, einsum, nan_to_num, sinc, zeros_like, cross, arctan2, tanh, cos, sin, real, abs
from .constants import *
import numpy as np

def abphase2d(mx, my, mz, Lx=1e-6, Ly=1e-6, p=array([0,0,1]), t=60e-9):
	Mx = fftshift(fft2(fftshift(mx)))
	My = fftshift(fft2(fftshift(my)))
	Mz = fftshift(fft2(fftshift(mz)))
	M = array([Mx, My, Mz])
	xres, yres = int(M.shape[1]/2), int(M.shape[2]/2)
	SI = arange(-xres,xres)
	SJ = arange(-yres,yres)
	si, sj = meshgrid(SI, SI)
	s = array([si/Lx,sj/Ly,0*si])
	s_mag = sqrt(einsum('ijk,ijk->jk',s,s))
	sig = nan_to_num(s/s_mag, nan=0,posinf=0,neginf=0)
	Gts = 1/(einsum('i,ijk->jk',p,sig)**2 + p[2]**2) \
			* sinc(t*einsum('i,ijk->jk',p,sig)/p[2])
	weights = zeros_like(M[0], dtype=complex)

	pxpxM = cross(p, cross(p, M, axisb=0))
	sigxz = cross(sig, array([0,0,1]), axisa=0)
	d = einsum('ijk,ijk->ij', sigxz, pxpxM)
	weights = 1j * nan_to_num(t/s_mag, nan=0,neginf=0,posinf=0) * Gts * d

	phi = 2*e/hbar/c * fftshift(ifft2(fftshift(weights)))
	return(real(phi))

def propagate(x, y, cphase, defocus=0, wavelength=1.97e-12, focal_length=1):
	dx = x[1,1]-x[0,0]
	dy = y[1,1]-y[0,0]
	U = fftshift(fftfreq(x.shape[0],dx))
	V = fftshift(fftfreq(x.shape[1],dy))
	qx, qy = meshgrid(U, V)
	psi_0 = exp(1j*cphase)
	psi_q = fftshift(fft2(fftshift(psi_0)))
	psi_out = fftshift(ifft2(fftshift(psi_q * T(qx, qy, defocus, wavelength))))
	return(psi_out)

def T(qx, qy, defocus, wavelength):
	out = aperture(qx, qy) * exp(-1j*chi(qx, qy, defocus, wavelength))
	return(out)

def chi(qx, qy, defocus, wavelength):
	return(pi*wavelength*defocus*(qx**2 + qy**2))

def aperture(qx, qy):
	out = np.zeros_like(qx)
	out[sqrt(qx**2 + qy**2) < np.max(sqrt(qx**2 + qy**2))] = 1
	return(out)

def jchessmodel(x, y, z=0, **kwargs):
	p = {   'aa':5, 'ba':5, 'ca':0,
			'ak':5e7, 'bk':-5e1, 'ck':0,
			'bg':5e7, 'cg':pi/2, 'n': 1}
	for key in kwargs.keys():
		if not key in p.keys(): return("Error: {:} is not a kwarg.".format(key))
	p.update(kwargs)

	r, phi = sqrt(x**2+y**2), arctan2(y,x)

	alpha_z = p['aa'] * exp(-p['ba'] * z**2) + p['ca']
	k_z = p['ak'] * exp(-p['bk'] * z**2) + p['ck']
	gamma_z = pi / 2 * tanh(p['bg'] * z) + p['cg']
	Theta_rz = 2 * arctan2((k_z * r)**alpha_z, 1)

	mx = cos(p['n']*phi%(2*pi)-gamma_z) * sin(Theta_rz)
	my = sin(p['n']*phi%(2*pi)-gamma_z) * sin(Theta_rz)
	mz = cos(Theta_rz)
	return(mx, my, mz)
