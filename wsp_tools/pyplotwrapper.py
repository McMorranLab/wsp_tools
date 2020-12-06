"""A wrapper for matplotlib.pyplot containing common plotting routines.
"""
from . import plt, np, rgba

class singleAx():
	"""An extension of the matplotlib.axes.Axes class.

	This class adds macros for 2d plotting that I commonly use. In particular,
	it's easy to select only a window of your data to show, to add x-y axes,
	and to show the rgba version of a complex 2d array.

	Typical usage:

	```python
	x = np.linspace(-10,10,xres)
	y = np.linspace(-10,10,yres)
	data = np.cos(x+y)
	window = [-3,7,1,4]

	fig, ax = plt.subplots()
	myax = wsp_tools.singleAx(ax)
	ax.setAxes(x, y, window)
	ax.set_xytitle('x','y','title')
	plt.show()
	```

	More commonly, this class is returned by ```wsp_tools.pyplotwrapper.subplots```.

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
		ax.setXY(x_axis, y_axis, window)
		x_windowed, y_windowed, data_windowed = ax.prePlot(data)
		ax.ax.SomeOtherMatplotlibPlottingRoutine(x_windowed, y_windowed, data_windowed)
		plt.show()
		```

		Returns:

		1. xout - numpy.1darray: windowed x-axis
		2. yout - numpy.1darray: windowed y-axis
		3. dout - numpy.2darray: windowed data

		Additionally, it sets the ax element's extent.
		"""
		if not self.hasAxes:
			self.x = np.linspace(0,100,data.shape[1])
			self.y = np.linspace(0,100,data.shape[0])
		if not self.hasWindow:
			self.xmin, self.xmax = 0, 100
			self.ymin, self.ymax = 0, 100
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

		x, y should be numpy.1darray objects.

		Window should be wrt the values of x and y, not their
		indices.
		If window is omitted, the full data and x-y axes will be used.

		"""
		self.hasAxes = True
		self.x, self.y = x, y
		if not window is None:
			self.hasWindow = True
			self.xmin = window[0]
			self.xmax = window[1]
			self.ymin = window[2]
			self.ymax = window[3]

	def setWindow(self, window=(0,100,0,100)):
		"""Applies a window to the singleAx object.

		Can be applied before or after setAxes - whichever was applied last
		will be used.

		window = (xmin, xmax, ymin, ymax) wrt the values of the axes, not
		their indices.
		"""
		self.hasWindow = True
		self.xmin = window[0]
		self.xmax = window[1]
		self.ymin = window[2]
		self.ymax = window[3]

	def set_title(self, title='', **kwargs):
		"""Sets the title of the plot.

		Takes all the same args as matplotlib.axes.Axes.set_title
		"""
		self.ax.set_title(title, **kwargs)

	def set_xlabel(self, xlabel='', **kwargs):
		"""Sets the xlabel of the plot.

		Takes all the same args as matplotlib.axes.Axes.set_xlabel
		"""
		self.ax.set_xlabel(xlabel, **kwargs)

	def set_ylabel(self, ylabel='', **kwargs):
		"""Sets the ylabel of the plot.

		Takes all the same kwargs as matplotlib.axes.Axes.set_ylabel
		"""
		self.ax.set_ylabel(ylabel, **kwargs)

	def set_xytitle(self, xlabel='', ylabel='', title='', **kwargs):
		"""Set the xlabel, ylabel, and title at the same time.

		Sets all three even if not all are given.
		Takes all the same kwargs as matplotlib.axes.Axes.set_xlabel,
		matplotlib.axes.Axes.set_ylabel, matplotlib.axes.Axes.set_title.
		Whatever you input will be applied to all three.

		For individual control, use singleAx.set_xlabel, singleAx.set_ylabel,
		or singleAx.set_title.
		"""
		self.ax.set_xlabel(xlabel, **kwargs)
		self.ax.set_ylabel(ylabel, **kwargs)
		self.ax.set_title(title, **kwargs)

	def imshow(self, data, step=1, **kwargs):
		"""Imshows the data.

		If a window has been applied to the plot before this is called, it
		will be used.

		data should be a numpy.2darray object.
		step: data[::step,::step] will be shown.
		"""
		x, y, data = self.prePlot(data, step)
		self.ax.imshow(data, extent=self.extent, origin='lower', **kwargs)

	def quiver(self, data, step=1, **kwargs):
		"""Shows a quiver plot of complex data.

		data should be a complex numpy.2darray object.
			If real, the y component will just be zero everywhere.
		step: data[::step,::step] will be shown.
		"""
		xr, yr, data = self.prePlot(data, step)
		self.ax.quiver(xr,yr,np.real(data),np.imag(data),**kwargs)

	def rgba(self, data, step=1, alpha='intensity', cmap='uniform', **kwargs):
		"""Shows an rgba interpretation of complex data.

		Color represents complex angle, and brightness represents either
		amplitude or intensity.

		data should be a complex numpy.2darray object.
			If real, the color will be uniform and only brightness will vary.

		alpha can represent either 'intensity' (default) or 'amplitude'.
		"""

		x, y, data = self.prePlot(data)
		data = rgba(data,alpha=alpha,cmap=cmap)
		self.ax.imshow(data, extent=self.extent, origin='lower')

	def inset(self, window=None, color='white', **kwargs):
		"""Plots a square box with vertices defined by window.

		Window = (xmin, xmax, ymin, ymax).

		Default color is white. Takes all the same kwargs as
		matplotlib.pyplot.plot().
		"""
		self.ax.plot(np.linspace(window[0], window[1], 100),
						np.zeros(100) + window[2], color=color, **kwargs)
		self.ax.plot(np.linspace(window[0], window[1],100),
						np.zeros(100)+window[3], color=color, **kwargs)
		self.ax.plot(np.zeros(100) + window[0],
						np.linspace(window[2], window[3], 100),
						color=color, **kwargs)
		self.ax.plot(np.zeros(100) + window[1],
						np.linspace(window[2], window[3], 100),
						color=color, **kwargs)

def subplots(rc=11, sharex=False, sharey=False, squeeze=False,
			subplot_kw=None, gridspec_kw=None, **fig_kw):
		"""Creates a fig, ax instance but replaces ax with singleAx.

		Behaves identically to matplotlib.pyplot.subplots(), but replaces each
		matplotlib.axes.Axes object with a wsp_tools.pyplotwrapper.singleAx
		object.

		Each wsp_tools.pyplotwrapper.singleAx object in turn behaves just like a
		normal Axes object, but with added methods.
		"""
		fig, ax = plt.subplots(nrows=rc//10, ncols=rc%10,
							sharex=sharex, sharey=sharey, squeeze=squeeze,
							subplot_kw=subplot_kw, gridspec_kw=gridspec_kw,
							tight_layout=True, **fig_kw)
		for i in range(ax.shape[0]):
			for j in range(ax.shape[1]):
				ax[i][j] = singleAx(ax[i][j])
		return(fig, np.array(ax))
