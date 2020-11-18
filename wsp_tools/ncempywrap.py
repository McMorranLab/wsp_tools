import ncempy.io.dm as dm
import ncempy.io.ser as ser
from numpy import mean, std
import os, pathlib
from matplotlib.pyplot import imsave

def dm3topng(fname, fdir=None, outname=None, outdir='.'):
	if outname is None:
		outname = os.path.splitext(fname)[0]
		print(outname)
	if fdir is None:
		fdir = os.getcwd()
		print(fdir)
	file = dm.dmReader(os.path.join(fdir, fname))
	data = file['data']
	avg = mean(data)
	stdev = std(data)
	vmin = avg - 2*stdev
	vmax = avg + 2*stdev
	print(os.path.join(outdir,outname))
	imsave(arr=data, vmin=vmin, vmax=vmax, fname=os.path.join(outdir, outname+'.png'))

def sertopng(fname, fdir=None, outname=None, outdir='.'):
	if outname is None:
		outname = os.path.splitext(fname)[0]
		print(outname)
	if fdir is None:
		fdir = os.getcwd()
		print(fdir)
	file = ser.serReader(os.path.join(fdir, fname))
	data = file['data']
	avg = mean(data)
	stdev = std(data)
	vmin = avg - 2*stdev
	vmax = avg + 2*stdev
	print(os.path.join(outdir,outname))
	imsave(arr=data, vmin=vmin, vmax=vmax, fname=os.path.join(outdir, outname+'.png'))

##################################################
# %%
# testing
