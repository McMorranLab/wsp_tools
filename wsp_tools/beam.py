import numpy as np
from scipy.special import eval_genlaguerre, factorial
from .constants import *

def E(T_eV):
    """kinetic energy [eV] -> total energy [J]"""
    return(T_eV * e + m_e * c**2)

def p(T_eV):
    """kinetic energy [eV] -> momentum [N s]"""
    return(1/c * np.sqrt(E(T_eV)**2 - (m_e * c**2)**2))

def dB(T_eV):
    """kinetic energy [eV] -> deBroglie wavelength [m]"""
    return(h / p(T_eV))

def k(T_eV):
    """kinetic energy [eV] -> wavenumber [m^-1]"""
    return(p(T_eV) / hbar)

def omega(T_eV):
    """kinetic energy [eV] -> angular frequency [rad s^-1]"""
    return(E(T_eV) / hbar)

def v_g(T_eV):
    """kinetic energy [eV] -> group velocity [m s^-1]"""
    return(c * k(T_eV) / np.sqrt(m_e**2 * c**2 / hbar**2 + k(T_eV)**2))

def v_p(T_eV):
    """kinetic energy [eV] -> phase velocity [m s^-1]"""
    return(omega(T_eV)/k(T_eV))

def bessel(x, y, z = 0, l = 0, kz0 = k(3e5),
            kperp = 0.0005*k(3e5), dkperp = 0, N = 30):
    """
    Creates a bessel beam by adding plane waves. The spectrum is a circle in
    k-space (k_z, dkperp + kperp * cos(theta), kperp * sin(theta)).
    """
    mode = x*y*z*0j
    w_0 = np.sqrt((kz0**2 + kperp**2)*c**2 + m_e**2*c**4/hbar**2)
    for n in range(N):
        phi_n = 2*np.pi*n/N
        kxn, kyn, kzn = kperp*np.cos(phi_n) + dkperp, kperp*np.sin(phi_n), kz0
        w_n = np.sqrt(w_0**2 + dkperp**2 * c**2 + 2*kxn*dkperp*c**2)
        mode += np.exp(1j*(-w_n*t + kxn*x + kyn*y + kzn*z + l*phi_n))
    mode *= 1/np.sqrt(N)
    return(mode)

def besselPacket(t = 0, l = 0,
        kres = 2**7, kmin = -3*k(3e5), kmax = 3*k(3e5),
        kz0 = k(3e5), kperp = .5*k(3e5), dkperp = 0,
        sig = 0.05*k(3e5)):
    """
    Creates a bessel beam by Fourier transforming a Gaussian spectrum.
    The spectrum is a gaussian centered on a circle in k-space
    (k_z, dkperp + kperp * cos(theta), kperp * sin(theta)). """
    # K = np.linspace(kmin, kmax, kres)
    kx, ky, kz = np.ogrid[kmin:kmax:1j*kres,kmin:kmax:1j*kres,kmin:kmax:1j*kres]
    w = c * np.sqrt((kx**2 + ky**2 + kz**2) + c**2 * m_e**2 /hbar**2)
    # cylindrical
    kr = np.sqrt(kx**2 + ky**2)
    mode = kz*0j
    phi = np.arctan2(ky,kx)
    spec = np.exp( -( (kr-kperp)**2 + (kz-kz0)**2 ) / (2*sig**2)) \
            * np.exp(1j*phi*l-1j*w*t)
    mode = np.fft.ifftshift(np.fft.ifftn(np.fft.fftshift(spec)))
    return(mode)

def zR(k, w0):
    """Rayleigh range as a function of k, beam waist"""
    lam = 2*pi/k
    return(pi * w0**2 / lam)

def w(z,w0,k):
    """Beam waist as a function of z, beam waist at z=0, and k"""
    return(w0 * np.sqrt(1 + (z/zR(k,w0))**2))

np.seterr(divide='ignore')
def R(z, w0, k):
    """Radius of curvature as a function of z, w0, k"""
    return(z + np.divide(zR(k,w0)**2,z))

def LG(x, y, z = 0, l = 0, p = 0, w_0 = 2e-6 * m, lam = 1.97e-12 * m):
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
    return(mode)
