import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import gamma
from wsp_tools.constants import *
import wsp_tools.beam as beam

def B(z, B0=1, a=.05):
    """Returns the Glaser glockenfeld"""
    return( B0/(1+(z/a)**2) )

def dBdz(z, B0=1, a=.05):
    """Returns the first derivative of the Glaser glockenfeld"""
    return( -2 * z * a**2 * B0/(a**2+z**2)**2 )

def W(z, B0=1, a=.05):
    """Returns the Larmor frequency of a Glaser glockenfeld"""
    return(e/(2*m_e) * B(z, B0, a))

def dWdz(z, B0=1, a=.05):
    """The first derivative of the Larmor frequency of a Glaser glockenfeld"""
    return(e/(2*m_e) * dBdz(z, B0, a))

def f(z, y, B0=1, a=.05, k=beam.k(3e5)):
    """System of differential equations for four gauge-invariant quantities"""
    vals = np.array([-m_e * W(z,B0,a) * dWdz(z,B0,a) * y[1] + dWdz(z,B0,a)*y[2],
                2 * y[3] /hbar /k,
                m_e * dWdz(z,B0,a) * y[1] + 2 * m_e * W(z,B0,a) * y[3] /hbar /k,
                2 * m_e * (y[0] - W(z,B0,a) * y[2]) /hbar /k])
    return(vals)

def jac(z, y, B0=1, a=.05, k=beam.k(3e5)):
    """Jacobian for f()"""
    vals = np.array([[0, -m_e * W(z,B0,a) * dWdz(z,B0,a), dWdz(z,B0,a), 0],
            [0, 0, 0, 2 /hbar /k],
            [0, m_e * dWdz(z,B0,a), 0, 2 * m_e * W(z,B0,a) /hbar /k],
            [2 * m_e /hbar /k, 0, -2 * m_e * W(z,B0,a) /hbar /k, 0]])
    return(vals)

def r02f(z, w0, k,l):
    """Free space solution for <r^2> of LG mode."""
    return(1/2 * beam.w(z,w0,k)**2 * (1 + np.abs(l)))

def Tp0f(z, w0, k, l):
    """Free space solution for <T_{perp}> of LG mode."""
    R = beam.R(z,w0,k)
    w = beam.w(z,w0,k)
    return(hbar**2 / m_e * (1+np.abs(l))/w**2 * (1 + (k * w**2/2/R)**2))

def Lz0f(z,w0,k,l):
    """Free space solution for <L_z> of LG mode."""
    return(l * hbar + 0*z)

def Gp0f(z, w0, k, l):
    """Free space solution for <G_{perp}> of LG mode."""
    R = beam.R(z,w0,k)
    w = beam.w(z,w0,k)
    return(hbar * k * w**2/2/R * (1+np.abs(l)))

def solve_odes_LG(z0,zf,B0=0,a=.005,k=beam.k(3e5),w0=2e-6,l=0,method='RK45'):
    """Solves f(). Ignore message saying RK45 doesn't make use of the jacobian.
    That's true, but it's ok. """
    r02 = r02f(z0,w0,k,l)
    Tp0 = Tp0f(z0,w0,k,l) + m_e*W(z0,B0,a)**2*r02/2 + W(z0,B0,a)*l*hbar
    Lz0 = Lz0f(z0,w0,k,l) + m_e * W(z0,B0,a) * r02
    Gp0 = Gp0f(z0,w0,k,l)
    y0 = [Tp0, r02, Lz0, Gp0]
    r1 = solve_ivp(f, [z0, zf], y0, method=method,
                jac=jac, max_step = (zf-z0)/2048, args=(B0,a,k))
    return(r1)
