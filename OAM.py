#!/usr/bin/env python

"""
This module contains a variety of functions useful for simulating vortex modes.
It contains

1. Functions for generating coordinates
divide(a,b) - handles division by zero: 1/0 = np.infty
coords(kwargs*) - generates both real and reciprocal space coordinates
uvPlane(x,y,kwargs*) - generates reciprocal coordinates given real space coords
xy(kwargs*) - generates real space coordinates

2. Functions for generating modes
LG(x,y,kwargs*) - generates a general LG mode (can take z as float or array)
bessel(x,y,kwargs*) - generates a general bessel mode (can take z as float or array)
besselPacket(kwargs*) - generates a general bessel packet (centered at xyz=0)

3. Functions specific to the log-polar sorter
phiU(x,y,kwargs*) - generates default unwrapper phase
phiC(u,v,kwargs*) - generates default corrector phase
phaseElements(x,y,kwargs*) - generates both phase elements
sort(x,y,mode,kwargs*) - sorts a given mode

4. Functions to help plotting
bin1d(axis,func,dt,max_l) - bins a 1d array into l-spectrum
plotMode(x_,y_,mode,kwargs*) - plots the intensity then the phase of the 2d mode
"""

import numpy as np
from scipy.special import eval_genlaguerre, jv, hermite, factorial
from scipy import signal
import time
import os
import matplotlib.pyplot as plt
import st_Bessel as bessel

### Coord generation

def divide(a,b):
    b = np.array([b])
    out = a*b*1.0
    out[np.where(b != 0)] = a/b[np.where(b != 0)]
    out[np.where(b == 0)] = np.infty
    return(out)

def coords(
        xymin = -10e-6, xymax = 10e-6, res=1000,
        lam = 2e-12, f = 1):
    x,y = xy(xymin,xymin,xymax,xymax,res,res)
    u,v = uvPlane(x,y,lam,f)
    return(x,y,u,v)

def uvPlane(
        x, y, lam = 2e-12, f = 1):
    res = len(x)
    dx = (np.max(x)-np.min(x))/res
    U = lam*f*np.fft.fftshift(np.fft.fftfreq(res,dx))
    u,v = np.meshgrid(U,U)
    return(u,v)

def xy(
        xmin = -10e-6, ymin = -10e-6, xmax = 10e-6,
        ymax = 10e-6, xres = 1000, yres = 1000):
    X = np.linspace(xmin,xmax,xres); Y = np.linspace(ymin,ymax,yres)
    return(np.meshgrid(X,Y))

### Mode gen

def LG(
        x, y, z = 0, l = 0, p = 0, w_0 = 2E-6, lam = 2e-12):
    k = 2*np.pi/lam
    Clp = np.sqrt(2*factorial(p)/np.pi/factorial(p+np.abs(l)))
    z_R = np.pi*w_0**2/lam
    R = (z + divide(z_R**2,z))
    w = w_0*np.sqrt(1 + (z/z_R)**2)
    gouy = (2*p + np.abs(l) + 1) * np.nan_to_num(np.arctan2(z,z_R))
    r = np.sqrt(x**2+y**2); theta = np.arctan2(y,x)
    mode = Clp/w*(r*np.sqrt(2)/w)**np.abs(l) \
            * np.exp(-r**2/w**2) * eval_genlaguerre(p,np.abs(l),(2*r**2)/(w**2)) \
            * np.exp(-1j*k*r**2/2/R**2) * np.exp(-1j*l*theta) * np.exp(1j*gouy)
    mode /= np.max(np.abs(mode))
    return(mode)

def bessel(
            x, y, z = 0, t = 0, theta = 0, kx0 = 0,
            ky0 = 0, kz0 = 2*np.pi/2e-12, krad = 40e-8*2*np.pi/2e-12,
            N = 30, l = 0):
    m = 9.1e-31; c = 3e8; hbar = 1e-34
    mode = x*y*z*0j
    for n in range(N):
        phi_n = 2*np.pi*n/N
        kxn,kyn,kzn = krad*np.cos(phi_n) + kx0, krad*np.sin(phi_n)*np.cos(theta) + ky0, krad*np.sin(phi_n)*np.sin(theta) + kz0
        w_n = c * np.sqrt(kxn**2 + kyn**2 + kzn**2 + c**2/hbar**2*m**2)
        mode += np.exp(1j*(-w_n*t + kxn*x + kyn*y + kzn*z + l*phi_n))
    mode *= 1/np.sqrt(N)
    return(mode)

def besselPacket(
        kres = 2**7, kmin = -3*2*np.pi/2e-12, kmax = 3*2*np.pi/2e-12,
        krad = .5*2*np.pi/2e-12, t = 0, theta = 0, kx0 = 0, ky0 = 0,
        kz0 = 2*np.pi/2e-12, N = 30, l = 0, sig = .05*2*np.pi/2e-12):
    m = 9.1e-31; c = 3e8; hbar = 1e-34
    K = np.linspace(kmin,kmax,kres)
    k = np.array(np.meshgrid(K,K,K))
    w = c*np.sqrt(np.einsum('ijkl,ijkl->jkl',k,k) + c**2/hbar**2*m**2)
    mode = k[0]*0j
    for n in range(N):
        phi_n = 2*np.pi*n/N
        kappa_n = np.array([krad*np.cos(phi_n) + kx0, krad*np.sin(phi_n)*np.cos(theta) + ky0, krad*np.sin(phi_n)*np.sin(theta) + kz0])
        phi_k = np.exp(-((k[0]-kappa_n[0])**2+(k[1]-kappa_n[1])**2+(k[2]-kappa_n[2])**2)/(2*sig**2) + 1j*l*phi_n - 1j*t*w)
        phi_r = np.fft.ifftshift(np.fft.ifftn(np.fft.fftshift(phi_k)))
        mode += phi_r
    mode *= 1/np.sqrt(N)
    return(mode)

### Sorting simulations

def phiU(
        x, y, k = 2*np.pi/2e-12, f = 1, a = 50e-6, b = 10e-6):
    return(k*a/f*(y*np.arctan2(y,x) + x - np.log(np.sqrt(x**2+y**2)/b)))

def phiC(
        u, v, k = 2*np.pi/2e-12, f = 1, a = 50e-6, b = 10e-6):
    return(-k*a*b/f*np.exp(-u/a)*np.cos(v/a))

def phaseElements(
        x, y, a = 10e-6, b = 50e-6, lam = 2e-12, f = 1):
    k = 2*np.pi/lam
    u,v = uvPlane(x,y,lam,f)
    return(phiU(x,y,k,f,a,b),phiC(u,v,k,f,a,b))

def sort(
        x, y, mode, a = 10e-6, b = 50e-6, lam = 2e-12,
        f = 1, phiU = True, phiC = True):
    phi_U,phi_C = phaseElements(x,y,a,b,lam,f)
    if type(phiU) != bool:
        phi_U = phiU
    if type(phiC) != bool:
        phi_C = phiC
    mode1 = mode*np.exp(1j*phi_U)
    mode2 = np.fft.ifftshift(np.fft.fft2(np.fft.fftshift(mode1)))
    mode3 = mode2*np.exp(1j*phi_C)
    mode4 = np.fft.fftshift(np.fft.ifft2(np.fft.ifftshift(mode3)))
    return(mode1,mode2,mode3,mode4)

### Analysis methods

def bin1d(
        axis, func, dt, max_l):
    max_l = 4
    bins = np.arange(-max_l,max_l+1)
    bind = np.zeros_like(bins)
    for i in range(len(bins)):
        temp = np.sum(func[(axis > -dt/2 + bins[i]*dt) & (axis < dt/2 + bins[i]*dt)])
        bind[i] = temp
    return(bind)

def plotMode(
        x_, y_, mode, figsize = (12,5), window = 0,
        title = "", colorbar = False):
    x = 1e6*x_
    y = 1e6*y_
    plt.figure(figsize=figsize); plt.suptitle(title,fontsize=20)
    extent = [np.min(x),np.max(x),np.min(y),np.max(y)]
    if window == 0 or window >= len(mode)/2:
        plt.subplot(121); plt.title("Intensity")
        plt.xlabel("x (µm)"); plt.ylabel("y (µm)")
        plt.imshow(np.abs(mode)**2,cmap='binary_r',extent=extent)
        if colorbar: plt.colorbar()
        plt.subplot(122); plt.title("Phase")
        plt.xlabel("x (µm)"); plt.ylabel("y (µm)")
        plt.imshow(np.angle(mode),extent=extent)
        if colorbar: plt.colorbar()
    else:
        xymin = int(len(mode)/2 - window)
        xymax = int(len(mode)/2 + window)
        coords1,coords2 = x[xymin:xymax,xymin:xymax],y[xymin:xymax,xymin:xymax]
        extent = [np.min(coords1),np.max(coords1),np.min(coords2),np.max(coords2)]
        plt.subplot(121); plt.title("Intensity")
        plt.xlabel("x (µm)"); plt.ylabel("y (µm)")
        plt.imshow(np.abs(mode[xymin:xymax,xymin:xymax])**2,cmap='binary_r',extent=extent)
        if colorbar: plt.colorbar()
        plt.subplot(122); plt.title("Phase")
        plt.xlabel("x (µm)"); plt.ylabel("y (µm)")
        plt.imshow(np.angle(mode[xymin:xymax,xymin:xymax]),extent=extent)
        if colorbar: plt.colorbar()
    plt.show()
