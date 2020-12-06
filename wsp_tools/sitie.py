"""Contains utilities for reconstructing phase and magnetization from Lorentz images.

The most common use case is to generate a lorentz object from a ```.dm3``` file.
Then one can analyze using high_pass(), sitie(), crop_pixel_counts(), etc.

Example:

```python
import ncempy.io.dm as dm

fdir = '/path/to/data'
fname = 'lorentzImage.dm3'
file = dm.dmReader(os.path.join(fdir, fname))

img = wsp_tools.lorentz(dm3file)
img.crop_pixel_counts()
img.high_pass()
img.blur()
img.sitie()

### plot img.Bx, img.By, img.phase, img.data, img.rawData, etc
```
"""
from . import dB, rgba, np, ndi, plt, os
from . import constants as _
np.seterr(divide='ignore')
import json

class lorentz:
	"""Class that contains sitie information about a lorentz image.

	Input:

	dm3file: a dictionary-like object with the following keys:

	* data: numpy.2darray() containing the electron counts
	* pixelSize: [scalar, scalar] containing the x and y pixel sizes
	* pixelUnit: [string, string] containing the unit of the pixel sizes
	"""
	def __init__(self, dm3file):
		self.rawData = dm3file['data']
		self.data = dm3file['data']
		self.pixelSize = dm3file['pixelSize'][0]
		self.pixelUnit = dm3file['pixelUnit'][0]
		self.x = np.arange(0,self.data.shape[1]) * self.pixelSize
		self.y = np.arange(0,self.data.shape[0]) * self.pixelSize
		self.metadata = {
							'pixelSize':float(dm3file['pixelSize'][0]),
							'pixelUnit':dm3file['pixelUnit'][0]
						}
		self.phase = None
		self.Bx, self.By = None, None

	def reset(self):
		"""Resets data to the rawData.
		"""
		self.data = self.rawData

	def sitie(self, defocus, wavelength=1.96e-12):
		"""Carries out phase and B-field reconstruction.

		Assigns phase, Bx, and By attributes.
		"""
		self.metadata.update({'defocus': defocus})
		dummy = sitie_RHS(self.data, defocus, wavelength)
		self.phase = np.real(inverse_laplacian(dummy, self.pixelSize))
		self.Bx, self.By = B_from_phase(self.phase, thickness=1)

	def crop_pixel_counts(self, sigma=10):
		"""Crops any pixel counts that are higher or lower than some std from avg.

		Sets those pixels to avg +/- sigma*std.
		"""
		self.metadata.update({'crop pixel sigma': sigma})
		self.data = crop_pixel_counts(self.data, sigma=sigma)

	def high_pass(self, sigma=20):
		"""Applies a high-pass filter to the image data.
		"""
		self.metadata.update({'high pass sigma': sigma})
		self.data = high_pass(self.data, sigma=sigma)

	def low_pass(self, sigma=50):
		"""Applies a low-pass filter to the image data.
		"""
		self.metadata.update({'low pass sigma': sigma})
		self.data = low_pass(self.data, sigma=sigma)

	def blur(self, sigma=5, mode='wrap', cval=0.0):
		"""Applies a Gaussian blur to the image data.
		"""
		self.metadata.update({'blur sigma': sigma})
		self.data = blur(self.data, sigma, mode, cval)

	def preview(self, window=(0,-1,0,-1)):
		"""Preview the image.

		window = (xmin, xmax, ymin, ymax)
		"""
		(xmin, xmax, ymin, ymax) = window
		if xmax == -1:
			xmax = self.data.shape[0] - 1
		if ymax == -1:
			ymax = self.data.shape[1] - 1
		fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(6,6*(ymax-ymin)/(xmax-xmin)), tight_layout=True)
		data = self.data[ymin:ymax, xmin:xmax]
		data = crop_pixel_counts(data, sigma=10)
		extent = [self.x[xmin],self.x[xmax],self.y[ymin],self.y[ymax]]
		ax.set_title("Intensity", fontsize=30)
		ax.set_xlabel("x (px)")
		ax.set_ylabel("y (px)")
		ax.imshow(data, origin="lower", extent=extent)
		plt.show()

	def saveMeta(self, outdir='', outname='metadata.json'):
		"""Save the metadata of the lorentz object to a file."""
		with open(os.path.join(outdir, outname), 'w') as fp:
			json.dump(self.metadata, fp, indent="")

############################## Pre-processing ####################
def blur(image, sigma=5, mode='wrap', cval=0.0):
	"""Applies a Gaussian filter to the image.
	"""
	return(ndi.gaussian_filter(image, sigma, mode=mode, cval=cval))

def crop_pixel_counts(image, sigma=10):
	"""Crops the pixel counts to avg +/- sigma*std.
	"""
	avg = np.mean(image)
	std = np.std(image)
	vmin = avg - sigma*std
	vmax = avg + sigma*std
	image[image < vmin] = vmin
	image[image > vmax] = vmax
	return(image)

def high_pass(image, sigma=50):
	"""Applies a high-pass filter to the image data."""
	X = np.linspace(-image.shape[0]/2, image.shape[0]/2, image.shape[0])
	Y = np.linspace(-image.shape[1]/2, image.shape[1]/2, image.shape[1])
	x, y = np.meshgrid(X, Y)
	g = 1-np.exp(-(x**2 + y**2)/2/sigma**2)
	fft = np.fft.ifftshift(np.fft.fft2(np.fft.fftshift(image)))
	out = np.fft.ifftshift(np.fft.ifft2(np.fft.fftshift(fft * g)))
	return(np.real(out))

def low_pass(image, sigma=50):
	"""Applies a low-pass filter to the image data."""
	X = np.linspace(-image.shape[0]/2, image.shape[0]/2, image.shape[0])
	Y = np.linspace(-image.shape[1]/2, image.shape[1]/2, image.shape[1])
	x, y = np.meshgrid(X, Y)
	g = np.exp(-(x**2 + y**2)/2/sigma**2)
	fft = np.fft.ifftshift(np.fft.fft2(np.fft.fftshift(image)))
	out = np.fft.ifftshift(np.fft.ifft2(np.fft.fftshift(fft * g)))
	return(np.real(out))

################################## SITIE #######################################
def B_from_phase(phase, thickness=1):
	"""Calculates the transverse B-field that would impart a specific phase."""
	dpdy, dpdx = np.gradient(phase)
	Bx = _.hbar/_.e/thickness * dpdy
	By = -_.hbar/_.e/thickness * dpdx
	return(Bx, By)

def SITIE(image, defocus, pixel_size, wavelength=1.97e-12):
	"""Reconstruct the phase from a defocussed image."""
	f = sitie_RHS(image, defocus, wavelength)
	phase = inverse_laplacian(f, pixel_size)
	return(np.real(phase))

def sitie_RHS(I, defocus, wavelength=dB(3e5)):
	return(2*_.pi/defocus * (1 - I/np.mean(I)))

def inverse_laplacian(f, pixel_size):
	QX = np.fft.fftfreq(f.shape[0], pixel_size)
	QY = np.fft.fftfreq(f.shape[1], pixel_size)
	qx,qy = np.meshgrid(QX,QY)

	f = np.fft.fft2(np.fft.fftshift(f))
	f = f/(qy**2 + qx**2+np.max(qx)/qx.shape[0]/100000)
	f[(qx==0)&(qy==0)] = 0
	f = np.fft.ifftshift(np.fft.ifft2(f))
	return(np.real(f))
