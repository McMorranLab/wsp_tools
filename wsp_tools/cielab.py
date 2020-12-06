"""Module to generate rgba data from numpy.2darrays.
"""
from . import np, plt
import numpy as np
from . import constants as _

def cielab_image(data, alpha = 'intensity'):
	"""Converts a 2d complex array to rgba data based on the cielab color space.

	Input:

	numpy.2darray(). dtype can be real or complex - if real, the output color
	will be constant. If complex, the color will reflect the complex angle.

	Returns:

	numpy.ndarray() with shape [data.shape[0], data.shape[1], 4]
	"""
	rgb_image_components = np.zeros(np.append(data.shape, 4),dtype=np.uint8)
	if alpha == 'intensity':
		value = np.absolute(data)**2 * 255
	else:
		value = np.absolute(data) * 255
	hue = (np.angle(data) + np.pi) / 2
	pi6 = np.pi/6
	def sin2(array, offset):
		return(np.sin(array-offset)**2)
	r = sin2(hue, .15*pi6) + 0.35 * sin2(hue, 3.15 * pi6)
	b = sin2(hue, 4.25 * pi6)
	g = .6*sin2(hue, 2*pi6) + 0.065 * sin2(hue, 5.05 * pi6) + 0.445*b + 0.33*r
	rgb_image_components[:,:,0] = (r * value).astype(np.uint8)
	rgb_image_components[:,:,1] = (g * value).astype(np.uint8)
	rgb_image_components[:,:,2] = (b * value).astype(np.uint8)
	rgb_image_components[:,:,3] = np.full(data.shape, fill_value=255, dtype=np.uint8)
	return(rgb_image_components)

def rgba(mode, cmap = 'uniform', alpha = 'intensity'):
	"""Converts a 2d complex array to rgba data.

	Input:

	dtype can be real or complex - if real, the output color
	will be constant. If complex, the color will reflect the complex angle.

	If cmap = 'uniform', the cielab color space will be used.

	Alpha can be either intensity or amplitude.

	Output:

	numpy.ndarray() with shape [data.shape[0], data.shape[1], 4].
	"""
	mode /= np.max(np.abs(mode))
	if cmap == 'uniform':
		out = cielab_image(mode, alpha)
		return(out)
	colormap = plt.cm.ScalarMappable(cmap=cmap)
	out = colormap.to_rgba(np.angle(mode))
	out[:,:,-1] = np.abs(mode)**2/np.max(np.abs(mode)**2)
	return(out)

def whatIsC():
	"""Used for testing the constants functionality of the module.
	"""
	return(_.c)
