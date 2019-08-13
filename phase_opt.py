#! usr/bin/env python3

import numpy as np
from scipy.optimize import minimize

### Parameter fitting

def error_func(params,simdphase,fct,x,y):
    return(np.sum((fct(params,x,y)-simdphase)**2))

def modelPhi(params,x,y):
    A,B = params
    r,t = np.sqrt(x**2+y**2),np.arctan2(y,x)
    r[np.where(r==0)] = 1e-12
    k = 2*np.pi/2E-12; foc = 1
    dummy = k*A/foc*(r*t*np.sin(t) + r*np.cos(t)*(1-np.log(r/B)))
    dummy -= np.min(dummy)
    return(dummy)

def makeIdeal(phase,x,y):
    a = 10e-6; b = 50e-6
    phase -= np.min(phase)
    fitParams = minimize(error_func,np.array([a,b]),(phase,modelPhi,x,y),method='TNC')
    A,B = fitParams.x
    idealPhase = modelPhi((A,B),x,y)
    return(idealPhase,A,B)

def error_func2(params,perc,x,y):
    n,m,o,offset = params
    return(np.sum((perc - (n*(x-offset)**2 + m*y**2 + o*x))**2))

def correct_astig(phase,x,y):
    xx = np.linspace(-1,1,len(x))
    x1,y1 = np.meshgrid(xx,xx)
    for i in range(4):
        idealPhase,a,b = makeIdeal(phase,x,y)
        percDiff = 100*(phase-idealPhase)/(np.max(idealPhase)-np.min(idealPhase))
        fitParams = minimize(error_func2,np.array([0,0,0,-.3]),(percDiff,x1,y1))
        n,m,o,offset = fitParams.x
        astig = (np.max(idealPhase)-np.min(idealPhase))/100*(n*(x1-offset)**2 + m*(y1)**2 + o*x1)
        phase = phase - astig
    return(idealPhase,phase,a,b,n,m,o,offset)
