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
# init_bcparams:  f,X,Y,Z,main_args   -> bcparams
# setBCs:          f,X,Y,Z,bcparams   -> f
# init_f:           X,Y,Z,main_args   -> f
# re_init:         f,X,Y,Z,bcparams   -> f, bcparams (both optional, typically only bcparams)
################################################################################

lam = 1.97e-12
k = 2*np.pi/lam
a = 10e-6
b = 50e-6
foc = 1
CE = 6.53e6

def init_bcparams(f,X,Y,Z,main_args):

    bcparams = ()
    return(bcparams)

def setBCs(g,X,Y,Z,bcparams):
    ### Von-Neumann conditions
    g[0,:,:] = g[1,:,:]; g[-1,:,:] = g[-2,:,:]
    g[:,0,:] = g[:,1,:]; g[:,-1,:] = g[:,-2,:]
    g[:,:,0] = g[:,:,1]; g[:,:,-1] = g[:,:,-2]
    ### Dirichlet

    return(g)

def init_f(X,Y,Z,main_args):
    f = np.random.random(X.shape)
    return(f)

def re_init(f,X,Y,Z,bcparams):

    return(f,bcparams)

################################################################################

def dblRes(V,x,y,z):
    # Initialize coordinates and new solution
    time1 = time.time()
    dx = .5*(x[1]-x[0])
    xres = len(x); yres = len(y); zres = len(z)
    xmin,xmax = np.min(x),np.max(x)
    ymin,ymax = np.min(y),np.max(y)
    zmin,zmax = np.min(z),np.max(z)
    xnew = np.linspace(xmin,xmax,2*xres-1,dtype='float32')
    ynew = np.linspace(ymin,ymax,2*yres-1,dtype='float32')
    znew = np.linspace(zmin,zmax,2*zres-1,dtype='float32')
    newSol = np.zeros((V.shape[0]*2,V.shape[1]*2,V.shape[2]*2),dtype='float32')
    newSol[::2,::2,::2] = V

    # First set of fill-ins
    avgSten = np.array([[[1,1],[1,1]],[[1,1],[1,1]]])
    newSol[1:-1:2,1:-1:2,1:-1:2] = 1/8*convolve(avgSten,newSol[::2,::2,::2],'valid')
    diag = dx*np.sqrt(2)

    sten1 = 1/diag**2*np.array([[[1,1],[1,1]],[[0,0],[0,0]]]) # newSol[::2,::2,::2]
    sten2 = 1/dx**2*np.array([[[0,0],[0,0]],[[0,0],[1,1]]]) # newSol[1::2,1::2,1::2]

    newSol[1:-1:2,1:-1:2,2::2] = 1/2*1/(1/dx**2+2/diag**2)\
        *(convolve(np.transpose(sten1,axes=(1,2,0)),newSol[::2,::2,::2],'valid')\
          +convolve(np.transpose(sten2,axes=(0,1,2)),newSol[1::2,1::2,1::2],'valid'))
    newSol[1:-1:2,2::2,1:-1:2] = 1/2*1/(1/dx**2+2/diag**2)\
        *(convolve(np.transpose(sten1,axes=(2,0,1)),newSol[::2,::2,::2],'valid')\
          +convolve(np.transpose(sten2,axes=(1,2,0)),newSol[1::2,1::2,1::2],'valid'))
    newSol[2::2,1:-1:2,1:-1:2] = 1/2*1/(1/dx**2+2/diag**2)\
        *(convolve(np.transpose(sten1,axes=(0,1,2)),newSol[::2,::2,::2],'valid')\
          +convolve(np.transpose(sten2,axes=(2,0,1)),newSol[1::2,1::2,1::2],'valid'))
    time2 = time.time()
    print(time2-time1)
    # Second set of fill-ins
    sten = np.array([[[1,1]]])
    conv1 = convolve(np.transpose(sten,axes=(2,0,1)),newSol[1::2,1:-1:2,2::2],'valid')
    conv2 = convolve(np.transpose(sten,axes=(1,2,0)),newSol[1:-1:2,1::2,2::2],'valid')
    conv3 = convolve(np.transpose(sten,axes=(0,1,2)),newSol[2::2,1:-1:2,1::2],'valid')
    conv4 = convolve(np.transpose(sten,axes=(1,2,0)),newSol[2::2,1::2,1:-1:2],'valid')
    conv5 = convolve(np.transpose(sten,axes=(2,0,1)),newSol[1::2,2::2,1:-1:2],'valid')
    conv6 = convolve(np.transpose(sten,axes=(0,1,2)),newSol[1:-1:2,2::2,1::2],'valid')
    conv7 = convolve(np.transpose(sten,axes=(2,0,1)),newSol[::2,2::2,2::2],'valid')
    conv8 = convolve(np.transpose(sten,axes=(1,2,0)),newSol[2::2,::2,2::2],'valid')
    conv9 = convolve(np.transpose(sten,axes=(0,1,2)),newSol[2::2,2::2,::2],'valid')

    newSol[2::2,2::2,1:-1:2] = 1/6*(conv5 + conv4 + conv9)
    newSol[2::2,1:-1:2,2::2] = 1/6*(conv1 + conv3 + conv8)
    newSol[1:-1:2,2::2,2::2] = 1/6*(conv2 + conv6 + conv7)

    # Cuts off the endpoints in each dimension
    # Automatically fills in remaining boundary points with VN conditions
    newSol = newSol[:-1,:-1,:-1]
    newSol[0,:,:] = newSol[1,:,:]; newSol[-1,:,:] = newSol[-2,:,:]
    newSol[:,0,:] = newSol[:,1,:]; newSol[:,-1,:] = newSol[:,-2,:]
    newSol[:,:,0] = newSol[:,:,1]; newSol[:,:,-1] = newSol[:,:,-2]

    time3 = time.time()
    print(time3-time2)
    print(newSol.shape)
    return(newSol,xnew,ynew,znew)

def update(g,X,Y,Z,bcparams):
    dx = X[1,1,1]-X[0,0,0]
    dy,dz = dx,dx
    stencil = np.array([[[0,0,0],[0,1/dx**2,0],[0,0,0]],[[0,1/dy**2,0],\
            [1/dz**2,0,1/dz**2],[0,1/dy**2,0]],[[0,0,0],[0,1/dx**2,0],[0,0,0]]])
    temp1 = 1/2*1/(1/dx**2+1/dy**2+1/dz**2)*convolve(stencil,g,'valid')
    error = temp1 - g[1:-1,1:-1,1:-1]
    g[1:-1,1:-1,1:-1] = temp1
    g = setBCs(g,X,Y,Z,bcparams)
    epsErr = np.mean(np.abs(error)**2)/np.mean(np.abs(g)**2)
    return(g,epsErr)

def main(filename,main_args,numDub=2,numsRuns=[2**7,2**6,2**5]):
    fname = filename + '.npz'
    # Make coordinates
    xmin,xmax = -1*10E-6, 3*10E-6
    ymin,ymax = -10E-6, 10E-6
    zmin,zmax = -4*10E-6, 4*10E-6
    xres = 2**6; yres = 2**5; zres = 2**7
    dx = (xmax-xmin)/xres; dy = (ymax-ymin)/yres; dz = (zmax-zmin)/zres
    x = np.linspace(xmin,xmax,xres,dtype='float32')
    y = np.linspace(ymin,ymax,yres,dtype='float32')
    z = np.linspace(zmin,zmax,zres,dtype='float32')
    X,Y,Z = np.meshgrid(x,y,z)

    f = init_f(X,Y,Z,main_args)
    bcparams = init_bcparams(f,X,Y,Z,main_args)

    # j gives resolution doubling; keyboard interrupt halts the simulation
        # and saves the approximation
    for j in range(numDub+1):
        epsErr = np.infty
        if j>0:
            f,x,y,z = dblRes(f,x,y,z)
            X,Y,Z = np.meshgrid(x,y,z)
            f,bcparams = re_init(f,X,Y,Z,bcparams)
        i = 0
        try:
            while epsErr>1E-8 and i<numsRuns[j]:
                i+= 1
                f,epsErr = update(f,X,Y,Z,bcparams)
                print(i," error :",epsErr)
        except KeyboardInterrupt:
            break
        # plt.imshow(f[:,:,int(len(z)/2)]); plt.colorbar(); plt.show()
    np.savez(fname,f,x,y,z)
    print('File saved: {0}'.format(fname))

if __name__ == "__main__":
    import sys
    import argparse
    parser = argparse.ArgumentParser(description='parse some arguments')
    parser.add_argument('-filename',default='dummy')
    parser.add_argument('-numDub',default=2)
    parser.add_argument('-numsRuns',default=[2**7,2**6,2**5],nargs='*')
    ############################################################################
    ### Add arguments for main() ###############################################
    ############################################################################
    args = parser.parse_args()
    main(args.filename,(),args.numDub,args.numsRuns)
