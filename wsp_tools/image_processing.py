from . import np

# %%
def high_pass(data, sigma = 7):
	Fdata = np.fft.fft2(data)
	Xfreqs = np.fft.fftfreq(data.shape[1], 1/data.shape[1])
	Yfreqs = np.fft.fftfreq(data.shape[0], 1/data.shape[0])
	xfreqs, yfreqs = np.meshgrid(Xfreqs, Yfreqs)

	g = 1 - np.exp(-(xfreqs**2 + yfreqs**2)/2/sigma**2)
	FFdata = np.fft.ifft2(g * Fdata)
	return(FFdata)

def low_pass(data, sigma = 100):
	Fdata = np.fft.fft2(data)
	Xfreqs = np.fft.fftfreq(data.shape[1], 1/data.shape[1])
	Yfreqs = np.fft.fftfreq(data.shape[0], 1/data.shape[0])
	xfreqs, yfreqs = np.meshgrid(Xfreqs, Yfreqs)

	g = np.exp(-(xfreqs**2 + yfreqs**2)/2/sigma**2)
	FFdata = np.fft.ifft2(g * Fdata)
	return(FFdata)

def clip_data(data, sigma = 5):
	avg = np.mean(data)
	stdev = np.std(data)
	vmin = avg - sigma*stdev
	vmax = avg + sigma*stdev
	data[data < vmin] = vmin
	data[data > vmax] = vmax
	return(data)

def shift_pos(data):
	return(data - np.min(data))

# %%
class ndap(np.ndarray):
	def __new__(cls, data, x=None, y=None):
		dummy = np.asarray(data).copy().view(cls)
		return(dummy)

	def __init__(self, data, x = None, y = None):
		self.x = x
		self.y = y
		self.isComplex = np.iscomplexobj(data)

	def high_pass(self, sigma = 7):
		if self.isComplex:
			self[:,:] = high_pass(self, sigma)
		else:
			self[:,:] = np.real(high_pass(self, sigma))
		return(self)

	def low_pass(self, sigma = 100):
		if self.isComplex:
			self[:,:] = low_pass(self, sigma)
		else:
			self[:,:] = np.real(low_pass(self, sigma))
		return(self)

	def clip_data(self, sigma = 5):
		self[:,:] = clip_data(self, sigma)
		return(self)

	def shift_pos(self):
		self[:,:] = shift_pos(self)
		return(self)
