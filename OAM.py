#!/usr/bin/env python

import numpy as np
from scipy.special import eval_genlaguerre, jv, hermite
from scipy import signal
import time
import os
import matplotlib.pyplot as plt

def LG(x, y, l = 0, p = 0, w_0 = 2E-6):
    r = np.sqrt(x**2+y**2); theta = np.arctan2(y,x)
    mode = r**np.abs(l) * np.exp(-(r/w_0)**2)\
            * np.exp(-1j*l*theta)\
            * eval_genlaguerre(p,np.abs(l),(2*r**2)/(w_0**2))
    mode /= np.max(np.abs(mode))
    return(mode)

def uvPlane(x,y,lam=2e-12,f=1):
    res = len(x)
    dx = (np.max(x)-np.min(x))/res
    U = lam*f*np.fft.fftshift(np.fft.fftfreq(res,dx))
    u,v = np.meshgrid(U,U)
    return(u,v)

def phaseElements(x,y,a=10e-6,b=50e-6,lam=2e-12,f=1):
    k = 2*np.pi/lam
    phiU = a*k/f*(y*np.arctan2(y,x) - x*np.log(np.sqrt(x**2+y**2)/b) + x)
    u,v = uvPlane(x,y,lam,f)
    phiC = -b*a*k/f*np.exp(-u/a)*np.cos(v/a)
    return(phiU,phiC)

def sort(x,y,mode,a=10e-6,b=50e-6,lam=2e-12,f=1,phiU=True,phiC=True):
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

def bin1d(axis,func,dt,max_l):
    max_l = 4
    bins = np.arange(-max_l,max_l+1)
    bind = np.zeros_like(bins)
    for i in range(len(bins)):
        temp = np.sum(func[(axis > -dt/2 + bins[i]*dt) & (axis < dt/2 + bins[i]*dt)])
        bind[i] = temp
    return(bind)

def plotMode(x_,y_,mode,figsize=(12,5),window=0,title="",colorbar=False):
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
