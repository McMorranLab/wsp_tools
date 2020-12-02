from .beam import dB
from .constants import *
from .cielab import rgba
import numpy as np
import scipy.ndimage as ndi
np.seterr(divide='ignore')
from matplotlib.pyplot import subplots, imshow, quiver, show, tight_layout, savefig
import os
import errno
from shutil import copy
import json

class lorentz:
	def __init__(self, dm3file):
		"""
		dictionary-like dm3file with required keys:
			- data
			- filename
			- pixelUnit
			- pixelSize
		"""
		self.rawData = dm3file['data']
		self.data = dm3file['data']
		self.fname = os.path.splitext(dm3file['filename'])[0]
		if dm3file['pixelUnit'] == '\u00b5m':
			self.pixelSize = dm3file['pixelSize'][0] * 1e-6
		elif dm3file['pixelUnit'] == 'nm':
			self.pixelSize = dm3file['pixelSize'][0] * 1e-9
		elif dm3file['pixelUnit'] == 'mm':
			self.pixelSize = dm3file['pixelSize'][0] * 1e-3
		else:
			self.pixelSize = dm3file['pixelSize'][0]
		self.metadata = {
							'pixelSize':float(dm3file['pixelSize'][0]),
							'pixelUnit':dm3file['pixelUnit'][0]
						}
		self.phase = None
		self.Bx, self.By = None, None

	def reset(self):
		self.data = self.rawData

	def sitie(self, defocus, wavelength=1.96e-12):
		self.metadata.update({'defocus': defocus})
		dummy = sitie_RHS(self.data, defocus, wavelength)
		self.phase = np.real(inverse_laplacian(dummy, self.pixelSize))
		self.Bx, self.By = B_from_phase(self.phase, thickness=1)

	def crop_pixel_counts(self, sigma=10):
		self.metadata.update({'crop pixel sigma': sigma})
		self.data = crop_pixel_values(self.data, sigma=sigma)

	def high_pass(self, sigma=20):
		self.metadata.update({'high pass sigma': sigma})
		self.data = high_pass(self.data, sigma=sigma)

	def low_pass(self, sigma=50):
		self.metadata.update({'low pass sigma': sigma})
		self.data = low_pass(self.data, sigma=sigma)

	def blur(self, sigma=5, mode='wrap', cval=0.0):
		self.metadata.update({'blur sigma': sigma})
		self.data = blur(self.data, sigma, mode, cval)

	def preview(self, window=((0,-1),(0,-1))):
		((xmin, xmax), (ymin, ymax)) = window
		if xmax == -1:
			xmax = self.data.shape[0]
		if ymax == -1:
			ymax = self.data.shape[1]
		fig, ax = subplots(nrows=1, ncols=1, figsize=(6,6*(ymax-ymin)/(xmax-xmin)))
		data = self.data[ymin:ymax, xmin:xmax]
		extent = [xmin,xmax,ymin,ymax]
		ax.set_title("Intensity", fontsize=30)
		ax.set_xlabel("x (px)")
		ax.set_ylabel("y (px)")
		ax.imshow(data, origin="lower", extent=extent)
		tight_layout()
		show()

	def show(self, window=((0,-1), (0, -1)), quiver_step=32):
		if self.phase is None:
			self.sitie(1e-3)
		((xmin, xmax), (ymin, ymax)) = window
		if xmax == -1:
			xmax = self.data.shape[0]
		if ymax == -1:
			ymax = self.data.shape[1]
		extent=[xmin,xmax,ymin,ymax]
		xrange = np.arange(xmin,xmax,quiver_step)
		yrange = np.arange(ymin,ymax,quiver_step)
		data = self.data[ymin:ymax, xmin:xmax]
		phase = self.phase[ymin:ymax, xmin:xmax]
		Bx, By = self.Bx[ymin:ymax, xmin:xmax], self.By[ymin:ymax, xmin:xmax]
		fig, ax = subplots(nrows=3, ncols=1, sharex=True, sharey=True, figsize=(6,3*6*(ymax-ymin)/(xmax-xmin)))
		(ax1, ax2, ax3) = ax
		ax1.set_title("Intensity", fontsize=30)
		ax1.set_ylabel("y (px)")
		ax1.imshow(data, origin='lower',extent=extent)
		ax2.set_title("Phase", fontsize=30)
		ax2.set_ylabel("y (px)")
		ax2.imshow(phase, origin='lower', extent=extent)
		ax3.set_title("B Field", fontsize=30)
		ax3.set_xlabel("x (px)")
		ax3.set_ylabel("y (px)")
		ax3.imshow(rgba(Bx+1j*By), origin='lower', extent=extent)
		ax3.quiver(xrange, yrange, Bx[::quiver_step,::quiver_step], By[::quiver_step,::quiver_step], color='white')
		tight_layout()
		show()

	def save(self, window=((0,-1), (0, -1)), quiver_step=32, outdir=None, outname=None,
			savedm3=False, metadata={}, titles=True, axes=True, separate=True, vertical=True):
		if self.phase is None:
			self.sitie(1e-3)
		((xmin, xmax), (ymin, ymax)) = window
		if xmax == -1:
			xmax = self.data.shape[0]
		if ymax == -1:
			ymax = self.data.shape[1]
		self.metadata.update({"xmin":xmin,"xmax":xmax,"ymin":ymin,"ymax":ymax})
		extent = [xmin,xmax,ymin,ymax]
		xrange = np.arange(xmin,xmax,quiver_step)
		yrange = np.arange(ymin,ymax,quiver_step)
		data = self.data[ymin:ymax, xmin:xmax]
		phase = self.phase[ymin:ymax, xmin:xmax]
		Bx, By = self.Bx[ymin:ymax, xmin:xmax], self.By[ymin:ymax, xmin:xmax]

		if outdir is None:
			outdir = os.getcwd()
		if outname is None:
			outdir = os.path.join(outdir, self.fname+"_x_{:}_{:}_y_{:}_{:}".format(xmin,xmax,ymin,ymax))
		else:
			outdir = os.path.join(outdir, outname)
		if not os.path.exists(outdir):
			try:
				os.makedirs(outdir, 0o700)
			except OSError as e:
				if e.errno != errno.EEXIST:
					raise

		if (window[0][0]!=0) or (window[0][1]!=-1) or (window[1][0]!=0) or (window[1][1]!=-1):
			fig, ax = subplots(nrows=1, ncols=1, figsize=(6,6*(ymax-ymin)/(xmax-xmin)))
			ax.imshow(self.data, origin='lower', extent=(0,self.data.shape[0],0,self.data.shape[1]))
			ax.plot(np.zeros(100)+xmin, np.arange(ymin,ymax), 'white')
			ax.plot(np.zeros(100)+xmax, np.arange(ymin,ymax), 'white')
			ax.plot(np.arange(xmin,xmax), np.zeros(100)+ymin, 'white')
			ax.plot(np.arange(xmin,xmax), np.zeros(100)+ymax, 'white')
			if axes:
				ax.set_xlabel("x (px)")
				ax.set_ylabel("y (px)")
			else:
				ax.set_xticks([])
				ax.set_yticks([])
			tight_layout()
			savefig(os.path.join(outdir, "full.png"))

		if separate:
			fig, ax = subplots(nrows=1, ncols=1, figsize=(6,6*(ymax-ymin)/(xmax-xmin)))
			if titles: ax.set_title("Intensity", fontsize=30)
			if axes:
				ax.set_xlabel("x (px)")
				ax.set_ylabel("y (px)")
			else:
				ax.set_xticks([])
				ax.set_yticks([])
			ax.imshow(data, origin='lower', extent=extent)
			tight_layout()
			savefig(os.path.join(outdir, "intensity.png"))

			fig, ax = subplots(nrows=1, ncols=1, figsize=(6,6*(ymax-ymin)/(xmax-xmin)))
			if titles: ax.set_title("Phase", fontsize=30)
			if axes:
				ax.set_xlabel("x (px)")
				ax.set_ylabel("y (px)")
			else:
				ax.set_xticks([])
				ax.set_yticks([])
			ax.imshow(phase, origin='lower', extent=extent)
			tight_layout()
			savefig(os.path.join(outdir, "phase.png"))

			fig, ax = subplots(nrows=1, ncols=1, figsize=(6,6*(ymax-ymin)/(xmax-xmin)))
			if titles: ax.set_title("B Field", fontsize=30)
			if axes:
				ax.set_xlabel("x (px)")
				ax.set_ylabel("y (px)")
			else:
				ax.set_xticks([])
				ax.set_yticks([])
			ax.imshow(rgba(Bx+1j*By), origin='lower', extent=extent)
			ax.quiver(xrange, yrange, Bx[::quiver_step,::quiver_step], By[::quiver_step,::quiver_step], color='white')
			tight_layout()
			savefig(os.path.join(outdir, "BField.png"))

		if not separate:
			if vertical:
				fig, (ax1,ax2,ax3) = subplots(nrows=3, ncols=1, figsize=(6,3*6*(ymax-ymin)/(xmax-xmin)))
			else:
				fig, (ax1,ax2,ax3) = subplots(nrows=1, ncols=3, figsize=(3*6*(ymax-ymin)/(xmax-xmin), 6))
			if titles: ax1.set_title("Intensity", fontsize=30)
			if axes:
				ax1.set_ylabel("y (px)")
			else:
				ax1.set_xticks([])
				ax1.set_yticks([])
			ax1.imshow(data, origin='lower', extent=extent)

			if titles: ax2.set_title("Phase", fontsize=30)
			if axes:
				ax2.set_ylabel("y (px)")
			else:
				ax2.set_yticks([])
				ax2.set_xticks([])
			ax2.imshow(phase, origin='lower', extent=extent)

			if titles: ax3.set_title("B Field", fontsize=30)
			if axes:
				ax3.set_xlabel("x (px)")
				ax3.set_ylabel("y (px)")
			else:
				ax3.set_xticks([])
				ax3.set_yticks([])
			ax3.imshow(rgba(Bx+1j*By), origin='lower', extent=extent)
			ax3.quiver(xrange, yrange, Bx[::quiver_step,::quiver_step], By[::quiver_step,::quiver_step], color='white')
			tight_layout()
			savefig(os.path.join(outdir, "combined.png"))

		if savedm3:
			copy(self.source, os.path.join(outdir, self.fname+'.dm3'))

		self.metadata.update(metadata)

		with open(os.path.join(outdir, 'metadata.json'), 'w') as fp:
			json.dump(self.metadata, fp, indent="")

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
	dpdy, dpdx = np.gradient(phase)
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
