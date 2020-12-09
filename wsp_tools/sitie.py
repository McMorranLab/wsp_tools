"""Contains utilities for reconstructing phase and magnetization from Lorentz images.

The most common use case is to generate a lorentz object from a ```.dm3``` file.
Then you can analyze using high_pass(), sitie(), clip_data(), etc.

Example:

```python
import ncempy.io.dm as dm
import wsp_tools as wt

fname = '/path/to/data.dm3'
dm3file = dm.dmReader(fname)

img = wt.lorentz(dm3file)
img.clip_data()
img.high_pass()
img.blur()
img.sitie()

### plot img.Bx, img.By, img.phase, img.data, img.rawData, etc
```
"""

# %%
from . import dB, rgba, np, plt, os
from .image_processing import *
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
		self.data = ndap(dm3file['data'])
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

	def fix_units(self, unit=None):
		if unit is None:
			if self.pixelUnit == 'nm':
				unit = 1e-9
			elif self.pixelUnit == 'mm':
				unit = 1e-3
			elif self.pixelUnit == 'µm':
				unit = 1e-6
			elif self.pixelUnit == 'm':
				unit = 1
		print("Setting unit from {:0.0e} meters to 1 meter".format(unit))
		self.pixelSize *= unit
		self.pixelUnit = 'm'
		self.x *= unit
		self.y *= unit
		self.metadata.update({
			'pixelSize': float(self.pixelSize),
			'pixelUnit': self.pixelUnit
		})
		return(None)

	def sitie(self, defocus = 1, wavelength=1.97e-12):
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
		self.metadata.update({'defocus': defocus, 'wavelength': wavelength})
		self.phase = ndap(SITIE(self.data, defocus, self.pixelSize, wavelength))
		self.Bx, self.By = [ndap(arr) for arr in B_from_phase(self.phase, thickness=1)]
		return(None)

	def preview(self, window=None):
		"""Preview the image.

		Note that unlike `pyplotwrapper`, window is in units of pixels.

		**Parameters**

		* **window** : _array-like, optional_ <br />
		Format is `window = (xmin, xmax, ymin, ymax)`. <br />
		Default is `window = (0, -1, 0, -1)`
		"""
		fig, ax = subplots(11)
		if not window is None:
			ax[0,0].setWindow(window)
		data = clip_data(self.data, sigma=10)
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

def SITIE(image, defocus, pixel_size = 1, wavelength=1.97e-12):
	"""Reconstruct the phase from a defocussed image.

	**Parameters**

	* **image** : _ndarray_ <br />
	the 2d image data. <br />

	* **defocus** : _number_ <br />

	* **pixel_size** : _number_ <br />

	* **wavelength** : _number, optional_ <br />
	Default is `wavelength = 1.97e-12` (the relativistic wavelength of a 300kV electron).

	"""
	rhs = sitie_RHS(image, defocus, wavelength)
	phase = inverse_laplacian(rhs, pixel_size)
	return(phase)

def sitie_RHS(I, defocus, wavelength=dB(3e5)):
	return(2 * _.pi / defocus / wavelength * (1 - I/np.mean(I)))

def inverse_laplacian(f, pixel_size):
	QX = np.fft.fftfreq(f.shape[1], pixel_size)
	QY = np.fft.fftfreq(f.shape[0], pixel_size)
	qx,qy = np.meshgrid(QX,QY)
	f = np.fft.fft2(f)
	f = np.nan_to_num(-f/(qy**2 + qx**2), posinf=0, neginf=0)
	f = np.fft.ifft2(f)
	return(np.real(f))
