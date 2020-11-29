from .beam import dB
from .constants import *
from .cielab import rgba
import numpy as np
import scipy.ndimage as ndi
np.seterr(divide='ignore')
from matplotlib.pyplot import subplots, imshow, quiver, show, tight_layout, savefig
import os
import ncempy.io.dm as dm
import errno
from shutil import copy
import json

class lorentz:
	def __init__(self, fdir, fname):
		dm3file = dm.dmReader(os.path.join(fdir, fname))
		self.rawData = dm3file['data']
		self.data = dm3file['data']
		self.fname = os.path.splitext(dm3file['filename'])[0]
		if dm3file['pixelUnit'] == 'µm':
			self.pixelSize = dm3file['pixelSize'][0] * 1e-6
		elif dm3file['pixelUnit'] == 'nm':
			self.pixelSize = dm3file['pixelSize'][0] * 1e-9
		elif dm3file['pixelUnit'] == 'mm':
			self.pixelSize = dm3file['pixelSize'][0] * 1e-3
		else:
			self.pixelSize = dm3file['pixelSize'][0]
		self.rawSITIE = SITIE(self.rawData, 1e-3, self.pixelSize)
		self.metadata = {'pixelSize':float(self.pixelSize)}
		self.sitie(1e-3)

	def sitie(self, defocus, wavelength=1.96e-12):
		self.metadata.update({'defocus': defocus})
		dummy = sitie_RHS(self.data, defocus, wavelength)
		self.phase = np.real(inverse_laplacian(dummy, self.pixelSize))
		self.Bx, self.By = B_from_phase(self.phase, thickness=1)

	def crop_pixel_counts(self, sigma=10):
		self.metadata.update({'crop pixel sigma': sigma})
		self.data = crop_pixel_values(self.data, sigma=sigma)

	def high_pass(self, sigma=50):
		self.metadata.update({'high pass sigma': sigma})
		self.data = high_pass(self.data, sigma=sigma)

	def low_pass(self, sigma=50):
		self.metadata.update({'low pass sigma': sigma})
		self.data = low_pass(self.data, sigma=sigma)

	def blur(self, sigma=5, mode='wrap', cval=0.0):
		self.metadata.update({'blur sigma': sigma})
		self.data = blur(self.data, sigma, mode, cval)

	def show(self, window=((0,-1), (0, -1)), quiver_step=32):
		((xmin, xmax), (ymin, ymax)) = window
		if xmax == -1:
			xmax = self.data.shape[0]
		if ymax == -1:
			ymax = self.data.shape[1]
		xrange = np.arange(xmin,xmax,quiver_step)
		yrange = np.arange(ymin,ymax,quiver_step)
		data = self.data[ymin:ymax, xmin:xmax]
		phase = self.phase[ymin:ymax, xmin:xmax]
		Bx, By = self.Bx[ymin:ymax, xmin:xmax], self.By[ymin:ymax, xmin:xmax]
		fig, ax = subplots(nrows=3, ncols=1, sharex=True, sharey=True, figsize=(6,18))
		(ax1, ax2, ax3) = ax
		ax1.set_title("Intensity - {:}".format(self.fname))
		ax1.set_ylabel("y (px)")
		ax1.imshow(data, origin='lower')
		ax2.set_title("Phase")
		ax2.set_ylabel("y (px)")
		ax2.imshow(phase, origin='lower')
		ax3.set_title("Phase")
		ax3.set_xlabel("x (px)")
		ax3.set_ylabel("y (px)")
		ax3.imshow(rgba(Bx+1j*By), origin='lower')
		ax3.quiver(xrange, yrange, By[::quiver_step,::quiver_step], Bx[::quiver_step,::quiver_step], color='white')
		tight_layout()
		show()

	def save(self, window=((0,-1), (0, -1)), quiver_step=32, outdir=None,
			savedm3=False, metadata={}):
		((xmin, xmax), (ymin, ymax)) = window
		if xmax == -1:
			xmax = self.data.shape[0]
		if ymax == -1:
			ymax = self.data.shape[1]
		xrange = np.arange(xmin,xmax,quiver_step)
		yrange = np.arange(ymin,ymax,quiver_step)
		data = self.data[ymin:ymax, xmin:xmax]
		phase = self.phase[ymin:ymax, xmin:xmax]
		Bx, By = self.Bx[ymin:ymax, xmin:xmax], self.By[ymin:ymax, xmin:xmax]
		if outdir is None:
			outdir = os.getcwd()
		outdir = os.path.join(outdir, self.fname+"_x_{:}_{:}_y_{:}_{:}".format(xmin,xmax,ymin,ymax))
		if not os.path.exists(outdir):
			try:
				os.makedirs(outdir, 0o700)
			except OSError as e:
				if e.errno != errno.EEXIST:
					raise
		fig, ax = subplots(nrows=1, ncols=1, figsize=(6,6))
		ax.set_title("Intensity - {:}".format(self.fname))
		ax.set_xlabel("x (px)")
		ax.set_ylabel("y (px)")
		ax.imshow(data, origin='lower')
		tight_layout()
		savefig(os.path.join(outdir, "intensity.png"))

		fig, ax = subplots(nrows=1, ncols=1, figsize=(6,6))
		ax.set_title("Phase - {:}".format(self.fname))
		ax.set_xlabel("x (px)")
		ax.set_ylabel("y (px)")
		ax.imshow(phase, origin='lower')
		tight_layout()
		savefig(os.path.join(outdir, "phase.png"))

		fig, ax = subplots(nrows=1, ncols=1, figsize=(6,6))
		ax.set_title("B Field - {:}".format(self.fname))
		ax.set_xlabel("x (px)")
		ax.set_ylabel("y (px)")
		ax.imshow(rgba(Bx+1j*By), origin='lower')
		ax.quiver(xrange, yrange, By[::quiver_step,::quiver_step], Bx[::quiver_step,::quiver_step], color='white')
		tight_layout()
		savefig(os.path.join(outdir, "BField.png"))

		if savedm3:
			copy(os.path.join(fdir, fname), os.path.join(outdir, fname))

		self.metadata.update(metadata)

		with open(os.path.join(outdir, 'metadata.json'), 'w') as fp:
			json.dump(self.metadata, fp)

############################## Pre-processing ####################
def blur(image, sigma=5, mode='wrap', cval=0.0):
	return(ndi.gaussian_filter(image, sigma, mode=mode, cval=cval))

def crop_pixel_values(image, sigma=10):
	avg = np.mean(image)
	std = np.std(image)
	vmin = avg - sigma*std
	vmax = avg + sigma*std
	image[image < vmin] = vmin
	image[image > vmax] = vmax
	return(image)

def high_pass(image, sigma=50):
	X = np.linspace(-image.shape[0]/2, image.shape[0]/2, image.shape[0])
	Y = np.linspace(-image.shape[1]/2, image.shape[1]/2, image.shape[1])
	x, y = np.meshgrid(X, Y)
	g = 1-np.exp(-(x**2 + y**2)/2/sigma**2)
	fft = np.fft.ifftshift(np.fft.fft2(np.fft.fftshift(image)))
	out = np.fft.ifftshift(np.fft.ifft2(np.fft.fftshift(fft * g)))
	return(np.real(out))

def low_pass(image, sigma=50):
	X = np.linspace(-image.shape[0]/2, image.shape[0]/2, image.shape[0])
	Y = np.linspace(-image.shape[1]/2, image.shape[1]/2, image.shape[1])
	x, y = np.meshgrid(X, Y)
	g = np.exp(-(x**2 + y**2)/2/sigma**2)
	fft = np.fft.ifftshift(np.fft.fft2(np.fft.fftshift(image)))
	out = np.fft.ifftshift(np.fft.ifft2(np.fft.fftshift(fft * g)))
	return(np.real(out))

################################## SITIE #######################################
def B_from_phase(phase, thickness=1):
	dpdx, dpdy = np.gradient(phase)
	Bx = hbar/e/thickness * dpdy
	By = -hbar/e/thickness * dpdx
	return(Bx, By)

def SITIE(image, defocus, pixel_size, wavelength=1.97e-12):
	f = sitie_RHS(image, defocus, wavelength)
	phase = inverse_laplacian(f, pixel_size)
	return(np.real(phase))

def sitie_RHS(I, defocus, wavelength=dB(3e5)):
	return(2*pi/defocus * (1 - I/np.mean(I)))

def inverse_laplacian(f, pixel_size):

	QX = np.fft.fftfreq(f.shape[0], pixel_size)
	QY = np.fft.fftfreq(f.shape[1], pixel_size)
	qx,qy = np.meshgrid(QX,QY)

	f = np.fft.fft2(np.fft.fftshift(f))
	f = f/(qy**2 + qx**2+np.max(qx)/qx.shape[0]/100000)
	f[(qx==0)&(qy==0)] = 0
	f = np.fft.ifftshift(np.fft.ifft2(f))
	return(np.real(f))
