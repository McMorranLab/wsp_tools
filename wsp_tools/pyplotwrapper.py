from . import plt, np, rgba

class singleFig():
	def __init__(self, ax, title='', xlabel='', ylabel=''):
		self.ax = ax
		self.hasAxes = False
		self.hasWindow = False

	def prePlot(self, data, step=1):
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
		self.hasAxes = True
		self.x, self.y = x, y
		if not window is None:
			self.hasWindow = True
			self.xmin = window[0]
			self.xmax = window[1]
			self.ymin = window[2]
			self.ymax = window[3]

	def setWindow(self, window=(0,100,0,100)):
		self.hasWindow = True
		self.xmin = window[0]
		self.xmax = window[1]
		self.ymin = window[2]
		self.ymax = window[3]

	def set_title(self, title=''):
		self.ax.set_title(title)

	def set_xlabel(self, xlabel=''):
		self.ax.set_xlabel(xlabel)

	def set_ylabel(self, ylabel=''):
		self.ax.set_ylabel(ylabel)

	def set_xytitle(self, xlabel='', ylabel='', title=''):
		self.ax.set_xlabel(xlabel)
		self.ax.set_ylabel(ylabel)
		self.ax.set_title(title)

	def imshow(self, data, step=1, **kwargs):
		x, y, data = self.prePlot(data, step)
		self.ax.imshow(data, extent=self.extent, origin='lower', **kwargs)

	def quiver(self, data, step=1, **kwargs):
		xr, yr, data = self.prePlot(data, step)
		self.ax.quiver(xr,yr,np.real(data),np.imag(data),**kwargs)

	def rgba(self, data, step=1, alpha='intensity', **kwargs):
		x, y, data = self.prePlot(data)
		data = rgba(data,alpha=alpha)
		self.ax.imshow(data, extent=self.extent)

	def inset(self, window=None, color='white'):
		self.ax.plot(np.linspace(window[0],window[1],100), np.zeros(100)+window[2], color=color)
		self.ax.plot(np.linspace(window[0],window[1],100), np.zeros(100)+window[3], color=color)
		self.ax.plot(np.zeros(100)+window[0], np.linspace(window[2],window[3],100), color=color)
		self.ax.plot(np.zeros(100)+window[1], np.linspace(window[2],window[3],100), color=color)

def subplots(
		nrows=1,ncols=1,sharex=False,sharey=False,squeeze=False,
		subplot_kw=None,gridspec_kw=None,**fig_kw):
		fig, ax = plt.subplots(nrows=nrows,ncols=ncols,
							sharex=sharex,sharey=sharey,squeeze=squeeze,
							subplot_kw=subplot_kw,gridspec_kw=gridspec_kw,
							tight_layout=True,**fig_kw)
		for i in range(ax.shape[0]):
			for j in range(ax.shape[1]):
				ax[i][j] = singleFig(ax[i][j])
		return(fig, np.array(ax))
