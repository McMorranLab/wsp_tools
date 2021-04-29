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

"""A wrapper for matplotlib.pyplot containing common plotting routines.

Example

```python
import numpy as np
import wsp_tools as wt

X = np.linspace(-10,10,100)
Y = np.linspace(-10,10,100)
x, y = np.meshgrid(X, Y)
data = x + 1j*y
window = (-4, 4, -6, 6)

fig, ax = wt.subplots(12)
ax1, ax2 = ax[0,0], ax[0,1]
ax1.setAxes(X, Y)
ax1.imshow(np.abs(data))
ax1.inset(window=window)
ax2.setAxes(X, Y, window=window)
ax2.rgba(data)
wt.plt.show()
```
"""

import matplotlib.pyplot as plt
import numpy as np
from .cielab import rgba, cielab_cmap

__all__ = ['singleAx','subplots']

class singleAx():
	"""An extension of the `matplotlib.axes.Axes` class.

	This class adds macros for 2d plotting that I commonly use. In particular,
	it's easy to select only a window of your data to show, to add x-y axes,
	to add an inset, and to show the rgba version of a complex 2d array.

	Typical usage:

	```python
	X = np.linspace(-10,10,xres)
	Y = np.linspace(-10,10,yres)
	x, y = np.meshgrid(X, Y)
	data = x+1j*y
	window = [-3,7,1,4]

	fig, ax = plt.subplots()
	myax = wsp_tools.singleAx(ax)
	ax.setAxes(x, y, window)
	ax.set_xytitle('x','y','title')
	ax.rgba(data)
	plt.show()
	```

	More commonly, this class is returned by ```wsp_tools.pyplotwrapper.subplots```.

	**Parameters**

	* **ax** : _matplotlib.axes.Axes_ <br />
	"""
	def __init__(self, ax, title='', xlabel='', ylabel=''):
		self.ax = ax
		self.hasAxes = False
		self.hasWindow = False
		self.set_xytitle(xlabel, ylabel, title)

	def prePlot(self, data, step=1):
		"""Utility function that applies the axes and window before plotting.

		If you want to use a plotting function from matplotlib, you can use this
		function to get the windowed data:

		```python
		fig, axis = plt.subplots()
		ax = singleAx(axis)
		ax.setWindow(window)
		x_windowed, y_windowed, data_windowed = ax.prePlot(data)
		ax.ax.SomeOtherMatplotlibPlottingRoutine(x_windowed, y_windowed, data_windowed)
		plt.show()
		```

		**Parameters** :

		* **data** : _complex ndarray_ <br />
		The data to plot. Must be 2-dimensional.

		* **step** : _int_ <br />
		data will be returned as `data[::step,::step]` - particularly useful for
		quiver plots. <br />
		Default is `step = 1`.

		**Returns**

		* **xout** : _ndarray_ <br />
		A 1darray with x-coordinates - either the array set by `setAxes()`,
		or an array of length `data.shape[1]` from 0 to 100.

		* **yout** : _ndarray_ <br />
		A 1darray with y-coordinates - either the array set by `setAxes()`,
		or an array of length `data.shape[0]` from 0 to 100.

		* **dout** : _ndarray_ <br />
		A 2darray with the data to be plotted. If you have set a window using
		either `setAxes()` or `setWindow()`, the data will be windowed.
		"""
		if not self.hasAxes and not self.hasWindow:
			self.x = np.linspace(0,100,data.shape[1])
			self.y = np.linspace(0,100,data.shape[0])
			self.argxmin, self.argxmax = 0, -1
			self.argymin, self.argymax = 0, -1
		elif self.hasAxes and not self.hasWindow:
			self.argxmin, self.argxmax = 0, -1
			self.argymin, self.argymax = 0, -1
		elif not self.hasAxes and self.hasWindow:
			self.x = np.linspace(0,100,data.shape[1])
			self.y = np.linspace(0,100,data.shape[0])
			self.argxmin = np.argmin(np.abs(self.x - self.xmin))
			self.argxmax = np.argmin(np.abs(self.x - self.xmax))
			self.argymin = np.argmin(np.abs(self.y - self.ymin))
			self.argymax = np.argmin(np.abs(self.y - self.ymax))
		elif self.hasAxes and self.hasWindow:
			self.argxmin = np.argmin(np.abs(self.x - self.xmin))
			self.argxmax = np.argmin(np.abs(self.x - self.xmax))
			self.argymin = np.argmin(np.abs(self.y - self.ymin))
			self.argymax = np.argmin(np.abs(self.y - self.ymax))

		xout = self.x[self.argxmin:self.argxmax:step]
		yout = self.y[self.argymin:self.argymax:step]
		dout = data[self.argymin:self.argymax:step, self.argxmin:self.argxmax:step]
		self.extent = [self.x[self.argxmin], self.x[self.argxmax], self.y[self.argymin], self.y[self.argymax]]
		return(xout, yout, dout)

	def setAxes(self, x, y, window=None):
		"""Sets the x and y axes of the singleAx object, and can apply a window.

		Note that this can be used before or after `setWindow()` - whichever
		was called last, will be used in plotting.

		**Parameters**

		* **x** : _ndarray_ <br />
		The x-coordinates. Should be 1-dimensional.

		* **y** : _ndarray_ <br />
		The y-coordinates. Should be 1-dimensional.

		* **window** : _array-like, optional_ <br />
		Format: `window = [xmin, xmax, ymin, ymax]`. Note that these are the x
		and y values, rather than their indices.

		**Returns**

		* **self** : _singleAx_
		"""
		self.hasAxes = True
		self.x, self.y = x, y
		if not window is None:
			self.hasWindow = True
			self.xmin = window[0]
			self.xmax = window[1]
			self.ymin = window[2]
			self.ymax = window[3]
		return(self)

	def setWindow(self, window=(0,100,0,100)):
		"""Applies a window to the singleAx object.

		Note that this can be used before or after `setAxes()` - whichever
		was called last, will be used in plotting.

		**Parameters**

		* **window** : _array-like, optional_ <br />
		Format: `window = [xmin, xmax, ymin, ymax]`. Note that these are the x
		and y values, rather than their indices.

		**Returns**

		* **self** : _singleAx_
		"""
		self.hasWindow = True
		self.xmin = window[0]
		self.xmax = window[1]
		self.ymin = window[2]
		self.ymax = window[3]
		return(self)

	def set_title(self, title='', **kwargs):
		"""Sets the title of the plot.

		**Parameters**

		* **title** : _string_ <br />
		The plot title.

		* ****kwargs** <br />
		All other kwargs are passed on to `matplotlib.axes.Axes.set_title`.

		**Returns**

		* **self** : _singleAx_
		"""
		self.ax.set_title(title, **kwargs)
		return(self)

	def set_xlabel(self, xlabel='', **kwargs):
		"""Sets the xlabel of the plot.

		**Parameters*

		* **xlabel** : _string_ <br />
		The xlabel.

		* ****kwargs** <br />
		All other kwargs are passed on to `matplotlib.axes.Axes.set_xlabel`.

		**Returns**

		* **self** : _singleAx_
		"""
		self.ax.set_xlabel(xlabel, **kwargs)
		return(self)

	def set_ylabel(self, ylabel='', **kwargs):
		"""Sets the ylabel of the plot.

		**Parameters**

		* **ylabel** : _string_ <br />
		The ylabel.

		* ****kwargs** <br />
		All other kwargs are passed on to `matplotlib.axes.Axes.set_ylabel`.

		**Returns**

		* **self** : _singleAx_
		"""
		self.ax.set_ylabel(ylabel, **kwargs)
		return(self)

	def set_xytitle(self, xlabel='', ylabel='', title='', **kwargs):
		"""Set the xlabel, ylabel, and title at the same time.

		Sets all three even if not all are given. Whatever you input will be applied to all three.

		For individual control, use `singleAx.set_xlabel`, `singleAx.set_ylabel`,
		or `singleAx.set_title`.

		**Parameters**

		* **ylabel** : _string_ <br />
		The ylabel.

		* **xlabel** : _string_ <br />
		The xlabel.

		* **title** : _string_ <br />
		The plot title.

		* ****kwargs** <br />
		All other kwargs are passed on
		to `matplotlib.axes.Axes.set_xlabel`, `matplotlib.axes.Axes.set_ylabel`,
		and `matplotlib.axes.Axes.set_title`.

		**Returns**

		* **self** : _singleAx_
		"""
		self.ax.set_xlabel(xlabel, **kwargs)
		self.ax.set_ylabel(ylabel, **kwargs)
		self.ax.set_title(title, **kwargs)
		return(self)

	def imshow(self, data, step=1, **kwargs):
		"""Imshows the (windowed) data. Sets origin to lower, sets the extent.

		**Parameters**

		* **data** : _ndarray_ <br />
		The data to be shown. Use the un-windowed data - the window will be
		applied automatically, if you set one.

		* **step** : _int_ <br />
		data will be returned as `data[::step,::step]` - particularly useful for
		quiver plots. <br />
		Default is `step = 1`.

		* ****kwargs** <br />
		All other kwargs are passed on to `matplotlib.axes.Axes.imshow`.

		**Returns**

		* **self** : _singleAx_
		"""
		x, y, data = self.prePlot(data, step)
		imshowargs = {'extent': self.extent, 'origin': 'lower'}
		imshowargs.update(kwargs)
		self.ax.imshow(data, **imshowargs)
		return(self)

	def quiver(self, data, step=1, rgba=False, **kwargs):
		"""Shows a quiver plot of complex data.

		**Parameters**

		* **data** : _ndarray_ <br />
		The data to be shown. Use the un-windowed data - the window will be
		applied automatically, if you set one.

		* **step** : _int_ <br />
		data will be returned as `data[::step,::step]` - particularly useful for
		quiver plots. <br />
		Default is `step = 1`.

		* **rgba** : _bool_ <br />
		If True, arrow color will correspond to the complex angle of the data. <br />
		Default is `rgba = False`.

		* ****kwargs** <br />
		All other kwargs are passed on to `matplotlib.axes.Axes.quiver`.

		**Returns**

		* **self** : _singleAx_
		"""
		xr, yr, data = self.prePlot(data, step)
		qargs = {'cmap' : cielab_cmap()}
		qargs.update(kwargs)
		self.ax.set_aspect('equal')
		if rgba:
			self.ax.set_facecolor('black')
			self.ax.quiver(xr, yr, np.real(data), np.imag(data), np.angle(data), **qargs)
		else:
			self.ax.quiver(xr,yr,np.real(data),np.imag(data),**kwargs)
		return(self)

	def rgba(self, data, step=1, brightness='intensity', alpha='uniform', cmap='uniform', **kwargs):
		"""Shows an rgba interpretation of complex data.

		**Parameters**

		* **data** : _complex ndarray_ <br />
		An array with the data to represent. Dtype may be complex or real - if real,
		the color will be uniform, and values will be represented by brightness.

		* **step** : _int_ <br />
		data will be returned as `data[::step,::step]` - particularly useful for
		quiver plots. <br />
		Default is `step = 1`.

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

		* ****kwargs** <br />
		All other kwargs are passed on to `matplotlib.axes.Axes.imshow`.

		**Returns**

		* **self** : _singleAx_
		"""

		x, y, data = self.prePlot(data)
		data = rgba(data,brightness=brightness,alpha=alpha,cmap=cmap)
		imshowargs = {'extent': self.extent, 'origin': 'lower'}
		imshowargs.update(kwargs)
		self.ax.imshow(data, **imshowargs)
		return(self)

	def inset(self, window, **kwargs):
		"""Plots a square box with vertices defined by window.

		Default color is white.

		**Parameters**

		* **window** : _array-like_ <br />
		Format: `window = [xmin, xmax, ymin, ymax]`. Note that these are the x
		and y values, rather than their indices.

		* ****kwargs** <br />
		All other kwargs are passed on to `matplotlib.axes.Axes.plot`.

		**Returns**

		* **self** : _singleAx_
		"""
		plotargs = {'color': 'white'}
		plotargs.update(kwargs)
		self.ax.plot(np.linspace(window[0], window[1], 100),
						np.zeros(100) + window[2], **plotargs)
		self.ax.plot(np.linspace(window[0], window[1],100),
						np.zeros(100)+window[3], **plotargs)
		self.ax.plot(np.zeros(100) + window[0],
						np.linspace(window[2], window[3], 100), **plotargs)
		self.ax.plot(np.zeros(100) + window[1],
						np.linspace(window[2], window[3], 100),
						**plotargs)
		return(self)

def subplots(rc=11, **kwargs):
		"""Creates a fig, ax instance but replaces ax with singleAx.

		Behaves almost identically to matplotlib.pyplot.subplots(), but replaces each
		`matplotlib.axes.Axes` object with a `wsp_tools.pyplotwrapper.singleAx`
		object.

		Each `wsp_tools.pyplotwrapper.singleAx` object in turn behaves just like a
		normal `Axes` object, but with added methods.

		**Parameters**

		* **rc** : _int_ <br />
		First digit - nrows. Second digit - ncols. <br />
		Default is `rc = 11`.

		* **squeeze** : _bool_ <br />
		If true, extra dimensions are squeezed out from the returned array of Axes. <br />
		Default is `squeeze = False`.

		* ****kwargs** <br />
		All other kwargs are passed on to `matplotlib.axes.Axes.subplots`.

		**Returns**

		* **fig** : _Figure_ <br />

		* **ax** : _singleAx_ or array of _singleAx_ objects <br />
		"""
		subplotsargs = {'tight_layout': True, 'squeeze': False}
		subplotsargs.update(kwargs)
		fig, ax = plt.subplots(nrows=rc//10, ncols=rc%10, **subplotsargs)
		for i in range(ax.shape[0]):
			for j in range(ax.shape[1]):
				ax[i][j] = singleAx(ax[i][j])
		return(fig, np.array(ax))
