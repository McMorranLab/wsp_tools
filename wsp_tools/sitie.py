from . import dB, rgba, np, ndi, plt, os
from .constants import *
np.seterr(divide='ignore')
import json

class lorentz:
	"""Class that contains sitie information about a lorentz image.

dictionary-like dm3file with required keys:
	- data
	- filename
	- pixelUnit
	- pixelSize
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
		self.data = self.rawData

	def sitie(self, defocus, wavelength=1.96e-12):
		self.metadata.update({'defocus': defocus})
		dummy = sitie_RHS(self.data, defocus, wavelength)
		self.phase = np.real(inverse_laplacian(dummy, self.pixelSize))
		self.Bx, self.By = B_from_phase(self.phase, thickness=1)

	def crop_pixel_counts(self, sigma=10):
		self.metadata.update({'crop pixel sigma': sigma})
		self.data = crop_pixel_counts(self.data, sigma=sigma)

	def high_pass(self, sigma=20):
		self.metadata.update({'high pass sigma': sigma})
		self.data = high_pass(self.data, sigma=sigma)

	def low_pass(self, sigma=50):
		self.metadata.update({'low pass sigma': sigma})
		self.data = low_pass(self.data, sigma=sigma)

	def blur(self, sigma=5, mode='wrap', cval=0.0):
		self.metadata.update({'blur sigma': sigma})
		self.data = blur(self.data, sigma, mode, cval)

	def preview(self, window=((0,-1),(0,-1))):
		((xmin, xmax), (ymin, ymax)) = window
		if xmax == -1:
			xmax = self.data.shape[0] - 1
		if ymax == -1:
			ymax = self.data.shape[1] - 1
		fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(6,6*(ymax-ymin)/(xmax-xmin)), tight_layout=True)
		data = self.data[ymin:ymax, xmin:xmax]
		data = crop_pixel_counts(data, sigma=10)
		extent = [self.x[xmin],self.x[xmax],self.y[ymin],self.y[ymax]]
		ax.set_title("Intensity", fontsize=30)
		ax.set_xlabel("x (px)")
		ax.set_ylabel("y (px)")
		ax.imshow(data, origin="lower", extent=extent)
		plt.show()

	def saveMeta(self, outdir=''):
		with open(os.path.join(outdir, 'metadata.json'), 'w') as fp:
			json.dump(self.metadata, fp, indent="")

############################## Pre-processing ####################
def blur(image, sigma=5, mode='wrap', cval=0.0):
	return(ndi.gaussian_filter(image, sigma, mode=mode, cval=cval))

def crop_pixel_counts(image, sigma=10):
	avg = np.mean(image)
	std = np.std(image)
	vmin = avg - sigma*std
	vmax = avg + sigma*std
	image[image < vmin] = vmin
	image[image > vmax] = vmax
	return(image)

def high_pass(image, sigma=50):
	X = np.linspace(-image.shape[0]/2, image.shape[0]/2, image.shape[0])
	Y = np.linspace(-image.shape[1]/2, image.shape[1]/2, image.shape[1])
	x, y = np.meshgrid(X, Y)
	g = 1-np.exp(-(x**2 + y**2)/2/sigma**2)
	fft = np.fft.ifftshift(np.fft.fft2(np.fft.fftshift(image)))
	out = np.fft.ifftshift(np.fft.ifft2(np.fft.fftshift(fft * g)))
	return(np.real(out))

def low_pass(image, sigma=50):
	X = np.linspace(-image.shape[0]/2, image.shape[0]/2, image.shape[0])
	Y = np.linspace(-image.shape[1]/2, image.shape[1]/2, image.shape[1])
	x, y = np.meshgrid(X, Y)
	g = np.exp(-(x**2 + y**2)/2/sigma**2)
	fft = np.fft.ifftshift(np.fft.fft2(np.fft.fftshift(image)))
	out = np.fft.ifftshift(np.fft.ifft2(np.fft.fftshift(fft * g)))
	return(np.real(out))

################################## SITIE #######################################
def B_from_phase(phase, thickness=1):
	dpdy, dpdx = np.gradient(phase)
	Bx = hbar/e/thickness * dpdy
	By = -hbar/e/thickness * dpdx
	return(Bx, By)

def SITIE(image, defocus, pixel_size, wavelength=1.97e-12):
	f = sitie_RHS(image, defocus, wavelength)
	phase = inverse_laplacian(f, pixel_size)
	return(np.real(phase))

def sitie_RHS(I, defocus, wavelength=dB(3e5)):
	return(2*pi/defocus * (1 - I/np.mean(I)))

def inverse_laplacian(f, pixel_size):

	QX = np.fft.fftfreq(f.shape[0], pixel_size)
	QY = np.fft.fftfreq(f.shape[1], pixel_size)
	qx,qy = np.meshgrid(QX,QY)

	f = np.fft.fft2(np.fft.fftshift(f))
	f = f/(qy**2 + qx**2+np.max(qx)/qx.shape[0]/100000)
	f[(qx==0)&(qy==0)] = 0
	f = np.fft.ifftshift(np.fft.ifft2(f))
	return(np.real(f))
