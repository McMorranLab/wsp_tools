"""Module to generate rgba data from scalar values.
"""
from . import np, plt
import numpy as np
from . import constants as _

def cielab_image(data, brightness = 'intensity', alpha = 'uniform'):
	"""Converts complex values to rgba data based on the CIELAB color space.

	The output color will represent the complex angle, and brightness may
	represent either intensity or amplitude.

	The CIELAB color space is intended to be perceptually uniform - none of the
	colors look brighter or darker than the others.

	**Parameters**

	* **data** : _ndarray_ <br />
	An array with the data to represent. Dtype may be complex or real - if real,
	the color will be uniform, and values will be represented by brightness.

	* **brightness** : _string, optional_ <br />
	Allowed values: `'intensity'`, `'amplitude'`, `'uniform'`. <br />
	Default is `brightness = 'intensity'`.

	* **alpha** : _string, optional_ <br />
	Allowed values: `'intensity'`, `'amplitude'`, `'uniform'`. Determines the alpha
	component of the rgba value. <br />
	Default is `alpha = 'uniform'`.

	**Returns**

	* **rgba_image_components** : _ndarray_ <br />
	The rgba components calculated from scalar values. If the input array has
	shape NxN, the output array will have shape NxNx4.
	"""
	data /= np.max(np.abs(data))
	rgba_image_components = np.zeros(np.append(data.shape, 4),dtype=np.uint8)
	if brightness == 'uniform':
		bvalue = 255
	elif brightness == 'intensity':
		bvalue = np.absolute(data)**2 * 255
	elif brightness == 'amplitude':
		bvalue = np.absolute(data) * 255
	if alpha == 'uniform':
		avalue = 255
	elif alpha == 'intensity':
		avalue = np.absolute(data)**2 * 255
	elif alpha == 'amplitude':
		avalue = np.absolute(data) * 255
	hue = (np.angle(data) + np.pi) / 2
	pi6 = np.pi/6
	def sin2(array, offset):
		return(np.sin(array-offset)**2)
	r = sin2(hue, .15*pi6) + 0.35 * sin2(hue, 3.15 * pi6)
	b = sin2(hue, 4.25 * pi6)
	g = .6*sin2(hue, 2*pi6) + 0.065 * sin2(hue, 5.05 * pi6) + 0.445*b + 0.33*r
	rgba_image_components[:,:,0] = (r * bvalue).astype(np.uint8)
	rgba_image_components[:,:,1] = (g * bvalue).astype(np.uint8)
	rgba_image_components[:,:,2] = (b * bvalue).astype(np.uint8)
	rgba_image_components[:,:,3] = np.full(data.shape, fill_value=avalue, dtype=np.uint8)
	return(rgba_image_components)

def rgba(mode, cmap = 'uniform', brightness = 'intensity', alpha = 'uniform'):
	"""Converts a 2d complex array to rgba data.

	**Parameters**

	* **mode** : _ndarray_ <br />
	An array with the data to represent. Dtype may be complex or real - if real,
	the color will be uniform, and values will be represented by brightness.

	* **cmap** : _string, optional_ <br />
	If `cmap = 'uniform'`, the CIELAB color space will be used. Otherwise, any
	pyplot ScalarMappable may be used. <br />
	Default is `cmap = 'uniform'`.

	* **brightness** : _string, optional_ <br />
	Allowed values: `'intensity'`, `'amplitude'`, `'uniform'`. <br />
	Default is `brightness = 'intensity'`.

	* **alpha** : _string, optional_ <br />
	Allowed values: `'intensity'`, `'amplitude'`, `'uniform'`. Determines the alpha
	component of the rgba value. <br />
	Default is `alpha = 'uniform'`.

	**Returns**

	* **rgba_image_components** : _ndarray_ <br />
	The rgba components calculated from scalar values. If the input array has
	shape NxN, the output array will have shape NxNx4.

	"""
	mode /= np.max(np.abs(mode))
	if cmap == 'uniform':
		out = cielab_image(mode, brightness, alpha)
		return(out)
	colormap = plt.cm.ScalarMappable(cmap=cmap)
	out = colormap.to_rgba(np.angle(mode))
	if alpha == 'intensity':
		out[...,-1] = np.abs(mode)**2
	elif alpha == 'amplitude':
		out[...,-1] = np.abs(mode)
	if brightness == 'intensity':
		out[...,0] *= np.abs(mode)**2
		out[...,1] *= np.abs(mode)**2
		out[...,2] *= np.abs(mode)**2
	elif brightness == 'amplitude':
		out[...,0] *= np.abs(mode)
		out[...,1] *= np.abs(mode)
		out[...,2] *= np.abs(mode)
	return(out)

def whatIsC():
	"""Used for testing the constants functionality of the module.
	"""
	return(_.c)
