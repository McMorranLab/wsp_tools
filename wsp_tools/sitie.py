"""Contains utilities for reconstructing phase and magnetization from Lorentz images.

The most common use case is to generate a lorentz object from a ```.dm3``` file.
Then you can analyze using high_pass(), sitie(), crop_pixel_counts(), etc.

Example:

```python
import ncempy.io.dm as dm
import wsp_tools as wt

fname = '/path/to/data.dm3'
dm3file = dm.dmReader(fname)

img = wt.lorentz(dm3file)
img.crop_pixel_counts()
img.high_pass()
img.blur()
img.sitie()

### plot img.Bx, img.By, img.phase, img.data, img.rawData, etc
```
"""
from . import dB, rgba, np, ndi, plt, os
from .pyplotwrapper import subplots
from . import constants as _
np.seterr(divide='ignore')
import json

class lorentz:
	"""Class that contains sitie information about a lorentz image.

	**Parameters**

	* **dm3file** : _dictionary-like_ <br />
	a dm3-like file with the following keys: <br />
		<ul>
		<li> **data** : _ndarray_ <br />
		An array carrying the electron counts. </li>
		<li> **pixelSize** : _tuple_ <br />
		(_number_, _number_) - the x and y pixel sizes. </li>
		<li> **pixelUnit** : _tuple_ <br />
		(_string_, _string_) - the x and y pixel units. </li>
		</ul>
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

		**Returns**

		* **self** : _lorentz_
		"""
		self.data = self.rawData
		return(self)

	def sitie(self, defocus, wavelength=1.97e-12):
		"""Carries out phase and B-field reconstruction.

		Assigns phase, Bx, and By attributes.

		**Parameters**

		* **defocus** : _number_ <br />
		The defocus at which the images were taken.

		* **wavelength** : _number, optional_ <br />
		The electron wavelength. <br />
		Default is `wavelength = 1.96e-12` (relativistic wavelength of a 300kV electron).

		**Returns**

		* **self** : _lorentz_
		"""
		self.metadata.update({'defocus': defocus})
		dummy = sitie_RHS(self.data, defocus, wavelength)
		self.phase = np.real(inverse_laplacian(dummy, self.pixelSize))
		self.Bx, self.By = B_from_phase(self.phase, thickness=1)
		return(self)

	def crop_pixel_counts(self, sigma=10):
		"""Crops any pixel counts that are higher or lower than some standard deviation from avg.

		All pixel counts are limited to (average) +/- sigma*(standard deviation).

		**Parameters**

		* **sigma** : _number, optional_ <br />
		Default is `sigma = 10`.

		**Returns**

		* **self** : _lorentz_
		"""
		self.metadata.update({'crop pixel sigma': sigma})
		self.data = crop_pixel_counts(self.data, sigma=sigma)
		return(self)

	def high_pass(self, sigma=20):
		"""Applies a high-pass filter to the image data.

		**Parameters**

		* **sigma** : _number, optional_ <br />
		Default is `sigma = 20`.

		**Returns**

		* **self** : _lorentz_
		"""
		self.metadata.update({'high pass sigma': sigma})
		self.data = high_pass(self.data, sigma=sigma)
		return(self)

	def low_pass(self, sigma=50):
		"""Applies a low-pass filter to the image data.

		**Parameters**

		* **sigma** : _number, optional_ <br />
		Default is `sigma = 50`.

		**Returns**

		* **self** : _lorentz_
		"""
		self.metadata.update({'low pass sigma': sigma})
		self.data = low_pass(self.data, sigma=sigma)
		return(self)

	def blur(self, sigma=5, mode='wrap', cval=0.0):
		"""Applies a Gaussian blur to the image data.

		**Parameters**

		* **sigma** : _number, optional_ <br />
		Default is `sigma = 5`.

		* **mode** : _string, optional_ <br />
		Gets passed to `scipy.ndimage.blur`. <br />
		Default is `mode = 'wrap'`.

		* **cval** : _number, optional_ <br />
		Gets passed to `scipy.ndimage.blur`. <br />
		Default is `cval = 0.0`.

		**Returns**

		* **self** : _lorentz_
		"""
		self.metadata.update({'blur sigma': sigma})
		self.data = blur(self.data, sigma, mode, cval)
		return(self)

	def preview(self, window=None):
		"""Preview the image.

		Note that unlike `pyplotwrapper`,

		**Parameters**

		* **window** : _array-like, optional_ <br />
		Format is `window = (xmin, xmax, ymin, ymax)`. <br />
		Default is `window = (0, -1, 0, -1)`
		"""
		fig, ax = subplots(11)
		if not window is None:
			ax[0,0].setWindow(window)
		data = crop_pixel_counts(self.data, sigma=10)
		ax[0,0].setAxes(self.x, self.y)
		ax[0,0].set_xlabel("x ({:})".format(self.pixelUnit))
		ax[0,0].set_ylabel("y ({:})".format(self.pixelUnit))
		ax[0,0].imshow(data)
		plt.show()

	def saveMeta(self, outdir='', outname='metadata.json'):
		"""Save the metadata of the lorentz object to a file.

		**Parameters**

		* **outdir** : _string, optional_ <br />
		The directory where you'd like to save the metadata. <br />
		Default is `outdir = ''`.

		* **outname** : _string, optional_ <br />
		The name of the metadata file. <br />
		Default is `outname = 'metadata.json'`.
		"""
		with open(os.path.join(outdir, outname), 'w') as fp:
			json.dump(self.metadata, fp, indent="")

############################## Pre-processing ####################
def blur(image, sigma=5, mode='wrap', cval=0.0):
	"""Applies a Gaussian filter to the image.

	**Parameters**

	* **sigma** : _number, optional_ <br />
	Default is `sigma = 5`.

	* **mode** : _string, optional_ <br />
	Gets passed to `scipy.ndimage.blur`. <br />
	Default is `mode = 'wrap'`.

	* **cval** : _number, optional_ <br />
	Gets passed to `scipy.ndimage.blur`. <br />
	Default is `cval = 0.0`.

	**Returns**

	* **None**
	"""
	return(ndi.gaussian_filter(image, sigma, mode=mode, cval=cval))

def crop_pixel_counts(image, sigma=10):
	"""Crops any pixel counts that are higher or lower than some standard deviation from avg.

	All pixel counts are limited to (average) +/- sigma*(standard deviation).

	**Parameters**

	* **sigma** : _number, optional_ <br />
	Default is `sigma = 10`.

	**Returns**

	* **None**
	"""
	avg = np.mean(image)
	std = np.std(image)
	vmin = avg - sigma*std
	vmax = avg + sigma*std
	image[image < vmin] = vmin
	image[image > vmax] = vmax
	return(image)

def high_pass(image, sigma=50):
	"""Applies a high-pass filter to the image data.

	**Parameters**

	* **sigma** : _number, optional_ <br />
	Default is `sigma = 20`.

	**Returns**

	* **None**
	"""
	X = np.linspace(-image.shape[0]/2, image.shape[0]/2, image.shape[0])
	Y = np.linspace(-image.shape[1]/2, image.shape[1]/2, image.shape[1])
	x, y = np.meshgrid(X, Y)
	g = 1-np.exp(-(x**2 + y**2)/2/sigma**2)
	fft = np.fft.ifftshift(np.fft.fft2(np.fft.fftshift(image)))
	out = np.fft.ifftshift(np.fft.ifft2(np.fft.fftshift(fft * g)))
	return(np.real(out))

def low_pass(image, sigma=50):
	"""Applies a low-pass filter to the image data.

	**Parameters**

	* **sigma** : _number, optional_ <br />
	Default is `sigma = 50`.

	**Returns**

	* **None**
	"""
	X = np.linspace(-image.shape[0]/2, image.shape[0]/2, image.shape[0])
	Y = np.linspace(-image.shape[1]/2, image.shape[1]/2, image.shape[1])
	x, y = np.meshgrid(X, Y)
	g = np.exp(-(x**2 + y**2)/2/sigma**2)
	fft = np.fft.ifftshift(np.fft.fft2(np.fft.fftshift(image)))
	out = np.fft.ifftshift(np.fft.ifft2(np.fft.fftshift(fft * g)))
	return(np.real(out))

################################## SITIE #######################################
def B_from_phase(phase, thickness=1):
	"""Reconstructs the B-field from the phase profile.

	**Parameters**

	* **phase** : _ndarray_ <br />
	a 2d array representing the electron's phase.

	* **thickness** : _number_ <br />
	the thickness of the sample. <br />
	Default is `thickness = 1`.

	**Returns**

	* **Bx** : _ndarray_ <br />
	The x-component of the magnetic field.

	* **By** : _ndarray_ <br />
	The y-component of the magnetic field.
	"""
	dpdy, dpdx = np.gradient(phase)
	Bx = _.hbar/_.e/thickness * dpdy
	By = -_.hbar/_.e/thickness * dpdx
	return(Bx, By)

def SITIE(image, defocus, pixel_size, wavelength=1.97e-12):
	"""Reconstruct the phase from a defocussed image.

	**Parameters**

	* **image** : _ndarray_ <br />
	the 2d image data. <br />

	* **defocus** : _number_ <br />

	* **pixel_size** : _number_ <br />

	* **wavelength** : _number, optional_ <br />
	Default is `wavelength = 1.97e-12` (the relativistic wavelength of a 300kV electron).

	"""
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
