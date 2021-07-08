# wsp-tools is TEM data analysis and simulation tools developed by WSP as a grad student in the McMorran Lab.
# Copyright (C) 2021  William S. Parker
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

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
img.sitie(defocus=1e-3)
img.phase.clip_data(sigma=5).high_pass().low_pass()
img.Bx.clip_data(sigma=5).high_pass().low_pass()
img.By.clip_data(sigma=5).high_pass().low_pass()

img.saveMeta(outdir='someDirectory') # Saves defocus, pixelSize, etc

### plot img.Bx, img.By, img.phase, img.data, img.rawData, etc
```
"""

# %%
import numpy as np
import matplotlib.pyplot as plt
import os
from deprecated import deprecated

from . import image_processing as ip
from .pyplotwrapper import subplots
from . import constants as _
from .ltem import phase_from_img, ind_from_phase
np.seterr(divide='ignore', invalid='ignore')
import json

__all__ = ['Lorentz','lorentz','B_from_phase','SITIE']

@deprecated(version='1.2.0', reason='This function has not found much use and will not be maintained')
def save_lorentz(img, fname=None, fdir=''):
	"""DEPRECATED. This function will no longer be maintained.

	Saves a `lorentz` object as a `.npz` archive.

	**Parameters**

	* **img** : _lorentz_ <br />
	The lorentz object to save.

	* **fname** : _string, optional_ <br />
	If not given, the output will be the filename in the lorentz object metadata, with `.npz` rather than `.dm3`.

	* **fdir** : _string, optional_ <br />
	The directory where the lorentz object will be saved. <br />
	Default is `fdir = ''`.

	**Returns**

	* **None**
	"""
	if fname is None:
		fname = os.path.splitext(img.metadata['filename'])[0]
	if not os.path.exists(fdir):
		if not fdir == '':
			os.makedirs(fdir)
	np.savez(os.path.join(fdir, fname), **img.__dict__)
	return(None)

@deprecated(version='1.2.0', reason='This function has not found much use and will not be maintained')
def load_lorentz(fname):
	"""DEPRECATED. This function will no longer be maintained.

	Loads a `lorentz` object that has been saved as a `.npz` via `save_lorentz()`.

	**Parameters**

	* **fname** : _string_ <br />

	**Returns**

	* **img** : _lorentz_ <br />
	The saved lorentz object, converted from a `.npz`.
	"""
	f = np.load(fname, allow_pickle=True)
	fm = f['metadata'].item()
	dm3 = {'data':f['data'], 'pixelSize':[fm['pixelSize'],fm['pixelSize']], 'pixelUnit':[fm['pixelUnit'],fm['pixelUnit']], 'filename':fm['filename']}
	img = lorentz(dm3)
	img.phase = ip.ndap(f['phase'])
	img.Bx = ip.ndap(f['Bx'])
	img.By = ip.ndap(f['By'])
	img.x = f['x']
	img.y = f['y']
	fm.update(img.metadata)
	img.metadata.update(fm)
	return(img)

def lorentz(file):
	"""Creates a `Lorentz` class instance for each image in a `.dm3` sequence.

	This function essentially acts as a wrapper for the `Lorentz` class, adding an extra handler
	for `.dm3` image sequences. When the `.dm3` file contains a sequence of images, it is split into
	individual `Lorentz` instances.

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
		<li> **filename** : _string_ <br /></li>
		</ul>

	**Returns**
	* **out** : _list, Lorentz_ <br />
	a list of `wsp-tools.sitie.Lorentz` instances, if the `.dm3` file is an image sequence.
	Otherwise, a single `wsp-tools.sitie.Lorentz` instance.
	"""
	if len(file['data'].shape) <= 2:
		return Lorentz(file)
	else:
		f1 = file.copy()
		out = []
		for dataset in file['data']:
			f1['data'] = dataset
			f1['pixelSize'] = [file['pixelSize'][-2], file['pixelSize'][-1]]
			f1['pixelUnit'] = [file['pixelUnit'][-2], file['pixelUnit'][-1]]
			f1['coords'] = [file['coords'][-2], file['coords'][-1]]
			out.append(Lorentz(f1))
		return(out)

class Lorentz:
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
		<li> **filename** : _string_ <br /></li>
		</ul>
	"""
	def __init__(self, dm3file):
		self.data = ip.ndap(dm3file['data'])
		self.dx = dm3file['pixelSize'][0]
		self.dy = dm3file['pixelSize'][1]
		self.xUnit = dm3file['pixelUnit'][0]
		self.yUnit = dm3file['pixelUnit'][1]
		self.x = np.arange(0,self.data.shape[1]) * self.dx
		self.y = np.arange(0,self.data.shape[0]) * self.dy
		self.metadata = {
							'dx':float(dm3file['pixelSize'][0]),
							'dy':float(dm3file['pixelSize'][1]),
							'xUnit':dm3file['pixelUnit'][0],
							'yUnit':dm3file['pixelUnit'][1],
							'filename':dm3file['filename']
						}
		self.phase = None
		self.Bx, self.By = None, None
		self.fix_units()

	def fix_units(self, xunit=None, yunit=None):
		"""Change the pixel units to meters.

		**Parameters**

		* **unit** : _number, optional_ <br />
		The scale to multiply values by (i.e., going from 'µm' to 'm', you would use `unit = 1e-6`). If none is given, `fix_units` will try to convert from `self.pixelUnit` to meters.

		**Returns**

		* **self** : _lorentz_
		"""
		if xunit is None:
			if self.xUnit == 'nm':
				xunit = 1e-9
			elif self.xUnit == 'mm':
				xunit = 1e-3
			elif self.xUnit == 'µm':
				xunit = 1e-6
			elif self.xUnit == 'm':
				xunit = 1
			else:
				xunit = 1
		if yunit is None:
			if self.yUnit == 'nm':
				yunit = 1e-9
			elif self.yUnit == 'mm':
				yunit = 1e-3
			elif self.yUnit == 'µm':
				yunit = 1e-6
			elif self.yUnit == 'm':
				yunit = 1
			else:
				yunit = 1
		self.dx *= xunit
		self.dy *= yunit
		self.xUnit = 'm'
		self.yUnit = 'm'
		self.x *= xunit
		self.y *= yunit
		self.metadata.update({
			'dx': float(self.dx),
			'dy': float(self.dy),
			'xUnit': self.xUnit,
			'yUnit': self.yUnit
		})
		return(None)

	def sitie(self, defocus = 1, thickness = 60e-9, wavelength=1.97e-12):
		"""Carries out phase and B-field reconstruction.

		Assigns phase, Bx, and By attributes.

		Updates metadata with the defocus and wavelength.

		**Parameters**

		* **defocus** : _number, optional_ <br />
		The defocus at which the images were taken. <br />
		Default is `defocus = 1`.

		* **wavelength** : _number, optional_ <br />
		The electron wavelength. <br />
		Default is `wavelength = 1.96e-12` (relativistic wavelength of a 300kV electron).

		**Returns**

		* **self** : _lorentz_
		"""
		self.metadata.update({'defocus': defocus, 'wavelength': wavelength, 'thickness': thickness})
		self.phase = ip.ndap(phase_from_img(self.data, defocus, self.dx, self.dy, wavelength))
		self.Bx, self.By = [ip.ndap(arr) for arr in ind_from_phase(self.phase, thickness)]
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
		data = ip.clip_data(self.data, sigma=10)
		ax[0,0].setAxes(self.x, self.y)
		ax[0,0].set_xlabel("x ({:})".format(self.dx))
		ax[0,0].set_ylabel("y ({:})".format(self.dy))
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
		if not os.path.exists(outdir):
			os.makedirs(outdir)
		with open(os.path.join(outdir, outname), 'w') as fp:
			json.dump(self.metadata, fp, indent="")

################################## SITIE #######################################
@deprecated(version='1.2.0', reason='You should now use wsp_tools.ltem.ind_from_phase. As of this version, B_from_phase just calls ind_from_phase.')
def B_from_phase(phase, thickness=1):
	"""DEPRECATED. You should instead use `wsp_tools.ltem.ind_from_phase`.

	Reconstructs the B-field from the phase profile.

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
	return(ind_from_phase(phase, thickness))

@deprecated(version='1.2.0', reason='You should now use wsp_tools.ltem.phase_from_img. As of this version, SITIE just calls phase_from_img.')
def SITIE(image, defocus, pixel_size = 1, wavelength=1.97e-12):
	"""DEPRECATED. You should instead use `wsp_tools.ltem.phase_from_img`.

	Reconstruct the phase from a defocussed image.

	**Parameters**

	* **image** : _ndarray_ <br />
	the 2d image data. <br />

	* **defocus** : _number_ <br />

	* **pixel_size** : _number, optional_ <br />
	Default is `pixel_size = 1`.

	* **wavelength** : _number, optional_ <br />
	Default is `wavelength = 1.97e-12` (the relativistic wavelength of a 300kV electron).

	**Returns**

	* **phase** : _ndarray_ <br />
	"""
	return(phase_from_img(img, defocus, pixel_size, pixel_size, wavelength))
