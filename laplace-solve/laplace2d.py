#! usr/bin/env python3

'''
Main arguments:
-filename: defaults to 'dummy.npz'
-numDubs: defaults to 2; number of times the resolution is doubled
-numsRuns: defaults to [2**7,2**6,2**5];
            the number of iterations to run at each resolution
'''

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import convolve
from scipy.optimize import minimize
import time

################################################################################
# Modify init_bcparams() to change how bcparams are initialized
# Modify setBCs() to change what boundary conditions are imposed
# Modify init_f() to change how f is initialized
# Modify re_init() to change what is reinitialized after every resolution doubling
# NOTE: Additional arguments can be easily added to main()

### Flow:
# init_bcparams:  f,X,Y,main_args   -> bcparams
# setBCs:          f,X,Y,bcparams   -> f
# init_f:           X,Y,main_args   -> f
# re_init:         f,X,Y,bcparams   -> f, bcparams (both optional, typically only bcparams)
################################################################################

def init_bcparams(f,X,Y,main_args):

    bcparams = ()
    return(bcparams)

def setBCs(g,X,Y,bcparams):
    ### Von-Neumann conditions
    g[0,:] = g[1,:]; g[-1,:] = g[-2,:]
    g[:,0] = g[:,1]; g[:,-1] = g[:,-2]
    ### Dirichlet
    ### BCs
    return(g)

def init_f(X,Y,main_args):
    f = np.random.random(X.shape)
    return(f)

def re_init(f,X,Y,bcparams):
    return(f,bcparams)

################################################################################

def dblRes(V,x,y):
    time1 = time.time()
    dx = x[1]-x[0]
    xres = len(x); yres = len(y)
    xmin,xmax = np.min(x),np.max(x)
    ymin,ymax = np.min(y),np.max(y)
    xnew = np.linspace(xmin,xmax,2*xres,dtype='float32')
    ynew = np.linspace(ymin,ymax,2*yres,dtype='float32')
    newSol = np.zeros((V.shape[0]*2,V.shape[1]*2),dtype='float32')
    newSol[::2,::2] = V
    sten1 = np.array([[1,1],[1,1]])
    newSol[1:-1:2,1:-1:2] = 1/4*convolve(sten1,newSol[::2,::2],'valid')
    sten2 = np.array([[1,1]])
    conv1 = convolve(np.transpose(sten2,axes=(1,0)),newSol[::2,2:-2:2],'valid')
    conv2 = convolve(np.transpose(sten2,axes=(0,1)),newSol[2:-2:2,::2],'valid')
    conv3 = convolve(np.transpose(sten2,axes=(0,1)),newSol[1:-1:2,1:-1:2],'valid')
    conv4 = convolve(np.transpose(sten2,axes=(1,0)),newSol[1:-1:2,1:-1:2],'valid')
    newSol[2:-2:2,1:-1:2] = 1/4*(conv2+conv4)
    newSol[1:-1:2,2:-2:2] = 1/4*(conv1+conv3)
    newSol = newSol[:-1,:-1]
    xnew = xnew[:-1]
    ynew = ynew[:-1]
    newSol[0,:] = newSol[1,:]; newSol[-1,:] = newSol[-2,:]
    newSol[:,0] = newSol[:,1]; newSol[:,-1] = newSol[:,-2]
    return(newSol,xnew,ynew)

def update(g,X,Y,bcparams):
    dx = X[1,1]-X[0,0]
    dy = dx
    stencil = np.array([[0,1,0],[1,0,1],[0,1,0]])
    temp1 = 1/4*convolve(stencil,g,'valid')
    error = temp1 - g[1:-1,1:-1]
    g[1:-1,1:-1] = temp1
    g = setBCs(g,X,Y,bcparams)
    epsErr = np.mean(np.abs(error)**2)/np.mean(np.abs(g)**2)
    return(g,epsErr)

def main(filename,main_args,numDub=2,numsRuns=[2**7,2**6,2**5]):
    fname = filename + '.npz'
    # Make coordinates
    xmin,xmax = -10, 10
    ymin,ymax = -10, 10
    xres = 2**6; yres = 2**6
    dx = (xmax-xmin)/xres; dy = (ymax-ymin)/yres
    x = np.linspace(xmin,xmax,xres,dtype='float32')
    y = np.linspace(ymin,ymax,yres,dtype='float32')
    X,Y = np.meshgrid(x,y)

    f = init_f(X,Y,main_args)
    bcparams = init_bcparams(f,X,Y,main_args)

    # j gives resolution doubling; keyboard interrupt halts the simulation
        # and saves the approximation
    for j in range(numDub+1):
        epsErr = np.infty
        if j>0:
            f,x,y = dblRes(f,x,y)
            X,Y = np.meshgrid(x,y)
            f,bcparams = re_init(f,X,Y,bcparams)
        i = 0
        try:
            while epsErr>1E-10 and i<numsRuns[j]:
                i+= 1
                f,epsErr = update(f,X,Y,bcparams)
                print(i," error :",epsErr)
        except KeyboardInterrupt:
            break
        plt.imshow(f); plt.colorbar(); plt.show()
    np.savez(fname,f,x,y)
    print('File saved: {0}'.format(fname))

if __name__ == "__main__":
    import sys
    import argparse
    parser = argparse.ArgumentParser(description='parse some arguments')
    parser.add_argument('-filename',default='dummy')
    parser.add_argument('-numDub',default=2)
    parser.add_argument('-numsRuns',default=[2**7,2**6,2**5],type=float,nargs='*')
    ############################################################################
    ### Add arguments for main() ###############################################
    ############################################################################
    args = parser.parse_args()
    print(args.numsRuns)
    main(args.filename,(),args.numDub,args.numsRuns)
