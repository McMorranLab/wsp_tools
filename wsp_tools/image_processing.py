from . import np

# %%
def high_pass(data, sigma = 7):
	"""Apply a high pass filter to a 2d-array.

	**Parameters**

	* **data** : _complex ndarray_ <br />

	* **sigma** : _number, optional_ <br />
	Standard deviation of the gaussian filter, measured in pixels. <br />
	Default is `sigma = 7`.

	**Returns**

	* **FFdata** : _complex ndarray_ <br />
	"""
	Fdata = np.fft.fft2(data)
	Xfreqs = np.fft.fftfreq(data.shape[1], 1/data.shape[1])
	Yfreqs = np.fft.fftfreq(data.shape[0], 1/data.shape[0])
	xfreqs, yfreqs = np.meshgrid(Xfreqs, Yfreqs)

	g = 1 - np.exp(-(xfreqs**2 + yfreqs**2)/2/sigma**2)
	FFdata = np.fft.ifft2(g * Fdata)
	return(FFdata)

def low_pass(data, sigma = 100):
	"""Apply a low pass filter to a 2d-array.

	**Parameters**

	* **data** : _complex ndarray_ <br />

	* **sigma** : _number, optional_ <br />
	Standard deviation of the gaussian filter, measured in pixels. <br />
	Default is `sigma = 100`.

	**Returns**

	* **FFdata** : _complex ndarray_ <br />
	"""
	Fdata = np.fft.fft2(data)
	Xfreqs = np.fft.fftfreq(data.shape[1], 1/data.shape[1])
	Yfreqs = np.fft.fftfreq(data.shape[0], 1/data.shape[0])
	xfreqs, yfreqs = np.meshgrid(Xfreqs, Yfreqs)

	g = np.exp(-(xfreqs**2 + yfreqs**2)/2/sigma**2)
	FFdata = np.fft.ifft2(g * Fdata)
	return(FFdata)

def clip_data(data, sigma = 5):
	"""Clip data to a certain number of standard deviations from average.

	* **data** : _complex ndarray_ <br />

	* **sigma** : _number, optional_ <br />
	Number of standard deviations from average to clip to. <br />
	Default is `sigma = 5`.

	**Returns**

	* **data** : _complex ndarray_ <br />
	"""
	avg = np.mean(data)
	stdev = np.std(data)
	vmin = avg - sigma*stdev
	vmax = avg + sigma*stdev
	data[data < vmin] = vmin
	data[data > vmax] = vmax
	return(data)

def shift_pos(data):
	"""Shift data to be all greater than zero.

	**Parameters**

	* **data** : _complex ndarray_ <br />

	**Returns**

	* **data** : _complex ndarray_
	"""
	return(data - np.min(data))

# %%
class ndap(np.ndarray):
	"""A class that adds all the image processing methods to np.ndarray.

	The purpose of this class is just so you can write `myarray.high_pass().low_pass()` instead of `myarray = high_pass(low_pass(myarray))`.

	**Parameters**

	* **data** : _complex ndarray_ <br />
	Any type of ndarray - the methods are defined with a 2d array in mind.
	"""
	def __new__(cls, data):
		dummy = np.asarray(data).copy().view(cls)
		return(dummy)

	def __init__(self, data):
		self.isComplex = np.iscomplexobj(data)

	def high_pass(self, sigma = 7):
		"""Apply a high pass filter to a 2d-array.

		**Parameters**

		* **sigma** : _number, optional_ <br />
		Standard deviation of the gaussian filter, measured in pixels. <br />
		Default is `sigma = 7`.

		**Returns**

		* **FFdata** : _ndap_ <br />
		"""
		if self.isComplex:
			self[:,:] = high_pass(self, sigma)
		else:
			self[:,:] = np.real(high_pass(self, sigma))
		return(self)

	def low_pass(self, sigma = 100):
		"""Apply a low pass filter to a 2d-array.

		**Parameters**

		* **sigma** : _number, optional_ <br />
		Standard deviation of the gaussian filter, measured in pixels. <br />
		Default is `sigma = 100`.

		**Returns**

		* **FFdata** : _ndap_ <br />
		"""
		if self.isComplex:
			self[:,:] = low_pass(self, sigma)
		else:
			self[:,:] = np.real(low_pass(self, sigma))
		return(self)

	def clip_data(self, sigma = 5):
		"""Clip data to a certain number of standard deviations from average.

		* **sigma** : _number, optional_ <br />
		Number of standard deviations from average to clip to. <br />
		Default is `sigma = 5`.

		**Returns**

		* **data** : _ndap_ <br />
		"""
		self[:,:] = clip_data(self, sigma)
		return(self)

	def shift_pos(self):
		"""Shift data to be all greater than zero.

		**Returns**

		* **data** : _ndap_
		"""
		self[:,:] = shift_pos(self)
		return(self)
