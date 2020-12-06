<a name="wsp_tools"></a>
# wsp\_tools

wsp_tools contains utilities for TEM data analysis and presentation.

Features:

* Single Image TIE
* Lorentz simulations
* spatial mode implementations
* a matplotlib.pyplot wrapper
* an implementation of the cielab colorspace
* a scipy.constants wrapper that allows unit scaling (i.e., using nanometers
instead of meters)

<a name="wsp_tools.sitie"></a>
# wsp\_tools.sitie

Contains utilities for reconstructing phase and magnetization from Lorentz images.

The most common use case is to generate a lorentz object from a ```.dm3``` file.
Then one can analyze using high_pass(), sitie(), crop_pixel_counts(), etc.

**Example**:

  
```python
import ncempy.io.dm as dm

fdir = '/path/to/data'
fname = 'lorentzImage.dm3'
file = dm.dmReader(os.path.join(fdir, fname))

img = wsp_tools.lorentz(dm3file)
img.crop_pixel_counts()
img.high_pass()
img.blur()
img.sitie()

### plot img.Bx, img.By, img.phase, img.data, img.rawData, etc
```

<a name="wsp_tools.sitie.lorentz"></a>
## lorentz Objects

```python
class lorentz()
```

Class that contains sitie information about a lorentz image.

Input:

* dm3file: a dictionary-like object with the following keys:
	* data: numpy.2darray() containing the electron counts
	* pixelSize: [scalar, scalar] containing the x and y pixel sizes
	* pixelUnit: [string, string] containing the unit of the pixel sizes

<a name="wsp_tools.sitie.lorentz.reset"></a>
#### reset

```python
 | reset()
```

Resets data to the rawData.

<a name="wsp_tools.sitie.lorentz.sitie"></a>
#### sitie

```python
 | sitie(defocus, wavelength=1.96e-12)
```

Carries out phase and B-field reconstruction.

Assigns phase, Bx, and By attributes.

<a name="wsp_tools.sitie.lorentz.crop_pixel_counts"></a>
#### crop\_pixel\_counts

```python
 | crop_pixel_counts(sigma=10)
```

Crops any pixel counts that are higher or lower than some std from avg.

Sets those pixels to avg +/- sigma*std.

<a name="wsp_tools.sitie.lorentz.high_pass"></a>
#### high\_pass

```python
 | high_pass(sigma=20)
```

Applies a high-pass filter to the image data.

<a name="wsp_tools.sitie.lorentz.low_pass"></a>
#### low\_pass

```python
 | low_pass(sigma=50)
```

Applies a low-pass filter to the image data.

<a name="wsp_tools.sitie.lorentz.blur"></a>
#### blur

```python
 | blur(sigma=5, mode='wrap', cval=0.0)
```

Applies a Gaussian blur to the image data.

<a name="wsp_tools.sitie.lorentz.preview"></a>
#### preview

```python
 | preview(window=(0,-1,0,-1))
```

Preview the image.

window = (xmin, xmax, ymin, ymax)

<a name="wsp_tools.sitie.lorentz.saveMeta"></a>
#### saveMeta

```python
 | saveMeta(outdir='', outname='metadata.json')
```

Save the metadata of the lorentz object to a file.

<a name="wsp_tools.sitie.blur"></a>
#### blur

```python
blur(image, sigma=5, mode='wrap', cval=0.0)
```

Applies a Gaussian filter to the image.

<a name="wsp_tools.sitie.crop_pixel_counts"></a>
#### crop\_pixel\_counts

```python
crop_pixel_counts(image, sigma=10)
```

Crops the pixel counts to avg +/- sigma*std.

<a name="wsp_tools.sitie.high_pass"></a>
#### high\_pass

```python
high_pass(image, sigma=50)
```

Applies a high-pass filter to the image data.

<a name="wsp_tools.sitie.low_pass"></a>
#### low\_pass

```python
low_pass(image, sigma=50)
```

Applies a low-pass filter to the image data.

<a name="wsp_tools.sitie.B_from_phase"></a>
#### B\_from\_phase

```python
B_from_phase(phase, thickness=1)
```

Calculates the transverse B-field that would impart a specific phase.

<a name="wsp_tools.sitie.SITIE"></a>
#### SITIE

```python
SITIE(image, defocus, pixel_size, wavelength=1.97e-12)
```

Reconstruct the phase from a defocussed image.

<a name="wsp_tools.pyplotwrapper"></a>
# wsp\_tools.pyplotwrapper

A wrapper for matplotlib.pyplot containing common plotting routines.

<a name="wsp_tools.pyplotwrapper.singleAx"></a>
## singleAx Objects

```python
class singleAx()
```

An extension of the matplotlib.axes.Axes class.

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

More commonly, this class is returned by
```wsp_tools.pyplotwrapper.subplots```.

<a name="wsp_tools.pyplotwrapper.singleAx.prePlot"></a>
#### prePlot

```python
 | prePlot(data, step=1)
```

Utility function that applies the axes and window before plotting.

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

**Returns**:

  
  1. xout - numpy.1darray: windowed x-axis
  2. yout - numpy.1darray: windowed y-axis
  3. dout - numpy.2darray: windowed data
  
  Additionally, it sets the ax element's extent.

<a name="wsp_tools.pyplotwrapper.singleAx.setAxes"></a>
#### setAxes

```python
 | setAxes(x, y, window=None)
```

Sets the x and y axes of the singleAx object, and can apply a window.

x, y should be numpy.1darray objects.

Window should be wrt the values of x and y, not their
indices.
If window is omitted, the full data and x-y axes will be used.

<a name="wsp_tools.pyplotwrapper.singleAx.setWindow"></a>
#### setWindow

```python
 | setWindow(window=(0,100,0,100))
```

Applies a window to the singleAx object.

Can be applied before or after setAxes - whichever was applied last
will be used.

window = (xmin, xmax, ymin, ymax) wrt the values of the axes, not
their indices.

<a name="wsp_tools.pyplotwrapper.singleAx.set_title"></a>
#### set\_title

```python
 | set_title(title='', **kwargs)
```

Sets the title of the plot.

Takes all the same args as matplotlib.axes.Axes.set_title

<a name="wsp_tools.pyplotwrapper.singleAx.set_xlabel"></a>
#### set\_xlabel

```python
 | set_xlabel(xlabel='', **kwargs)
```

Sets the xlabel of the plot.

Takes all the same args as matplotlib.axes.Axes.set_xlabel

<a name="wsp_tools.pyplotwrapper.singleAx.set_ylabel"></a>
#### set\_ylabel

```python
 | set_ylabel(ylabel='', **kwargs)
```

Sets the ylabel of the plot.

Takes all the same kwargs as matplotlib.axes.Axes.set_ylabel

<a name="wsp_tools.pyplotwrapper.singleAx.set_xytitle"></a>
#### set\_xytitle

```python
 | set_xytitle(xlabel='', ylabel='', title='', **kwargs)
```

Set the xlabel, ylabel, and title at the same time.

Sets all three even if not all are given.
Takes all the same kwargs as matplotlib.axes.Axes.set_xlabel,
matplotlib.axes.Axes.set_ylabel, matplotlib.axes.Axes.set_title.
Whatever you input will be applied to all three.

For individual control, use singleAx.set_xlabel, singleAx.set_ylabel,
or singleAx.set_title.

<a name="wsp_tools.pyplotwrapper.singleAx.imshow"></a>
#### imshow

```python
 | imshow(data, step=1, **kwargs)
```

Imshows the data.

If a window has been applied to the plot before this is called, it
will be used.

data should be a numpy.2darray object.
step: data[::step,::step] will be shown.

<a name="wsp_tools.pyplotwrapper.singleAx.quiver"></a>
#### quiver

```python
 | quiver(data, step=1, **kwargs)
```

Shows a quiver plot of complex data.

data should be a complex numpy.2darray object.
	If real, the y component will just be zero everywhere.
step: data[::step,::step] will be shown.

<a name="wsp_tools.pyplotwrapper.singleAx.rgba"></a>
#### rgba

```python
 | rgba(data, step=1, alpha='intensity', cmap='uniform', **kwargs)
```

Shows an rgba interpretation of complex data.

Color represents complex angle, and brightness represents either
amplitude or intensity.

data should be a complex numpy.2darray object.
	If real, the color will be uniform and only brightness will vary.

alpha can represent either 'intensity' (default) or 'amplitude'.

<a name="wsp_tools.pyplotwrapper.singleAx.inset"></a>
#### inset

```python
 | inset(window=None, color='white', **kwargs)
```

Plots a square box with vertices defined by window.

Window = (xmin, xmax, ymin, ymax).

Default color is white. Takes all the same kwargs as
matplotlib.pyplot.plot().

<a name="wsp_tools.pyplotwrapper.subplots"></a>
#### subplots

```python
subplots(rc=11, sharex=False, sharey=False, squeeze=False, subplot_kw=None, gridspec_kw=None, **fig_kw)
```

Creates a fig, ax instance but replaces ax with singleAx.

Behaves identically to matplotlib.pyplot.subplots(), but replaces each
matplotlib.axes.Axes object with a wsp_tools.pyplotwrapper.singleAx
object.

Each wsp_tools.pyplotwrapper.singleAx object in turn behaves just like a
normal Axes object, but with added methods.

<a name="wsp_tools.cielab"></a>
# wsp\_tools.cielab

Module to generate rgba data from numpy.2darrays.

<a name="wsp_tools.cielab.cielab_image"></a>
#### cielab\_image

```python
cielab_image(data, alpha='intensity')
```

Converts a 2d complex array to rgba data based on the cielab color space.

Input:

numpy.2darray(). dtype can be real or complex - if real, the output color
will be constant. If complex, the color will reflect the complex angle.

**Returns**:

  
  numpy.ndarray() with shape [data.shape[0], data.shape[1], 4]

<a name="wsp_tools.cielab.rgba"></a>
#### rgba

```python
rgba(mode, cmap='uniform', alpha='intensity')
```

Converts a 2d complex array to rgba data.

Input:

dtype can be real or complex - if real, the output color
will be constant. If complex, the color will reflect the complex angle.

If cmap = 'uniform', the cielab color space will be used.

Alpha can be either intensity or amplitude.

Output:

numpy.ndarray() with shape [data.shape[0], data.shape[1], 4].

<a name="wsp_tools.cielab.whatIsC"></a>
#### whatIsC

```python
whatIsC()
```

Used for testing the constants functionality of the module.

<a name="wsp_tools.constants"></a>
# wsp\_tools.constants

Wrapper for the scipy.constants module that allows unit scaling.

To see all available constants, print(wsp_tools.constants.__all__).

To set units:

import wsp_tools as wt
wt.setUnits(meter = 1e-3) # set km as base unit for length
from wsp_tools.constants import *

print(c) # outputs 299792.458

Note that all other modules should update automatically as well.

<a name="wsp_tools.constants.setUnits"></a>
#### setUnits

```python
setUnits(second=1, meter=1, kilogram=1, Amp=1, Kelvin=1, mole=1, candela=1)
```

Sets the units across the wsp_tools module.

i.e. setUnits(meter = 1000) sets the millimeter as the base unit for length.

<a name="wsp_tools.lorentzSim"></a>
# wsp\_tools.lorentzSim

<a name="wsp_tools.lorentzSim.abphase2d"></a>
#### abphase2d

```python
abphase2d(mx, my, mz, Lx=1e-6, Ly=1e-6, p=np.array([0,0,1]), t=60e-9)
```

Calculates the Aharonov-Bohm phase acquired by an electron through a 2d
magnetization.

Inputs:

The real space magnetization of the sample.

mx: numpy.2darray() - x component of magnetization

my: numpy.2darray() - y component of magnetization

mz: numpy.2darray() - z component of magnetization

Lx: scalar, the width of the sample.

Ly: scalar, the length of the sample.

p: numpy.1darray with length 3 - unit vector in the direction of electron motion.

t: the thickness of the sample.

**Returns**:

  
  phi = numpy.2darray() the phase acquired.

<a name="wsp_tools.lorentzSim.propagate"></a>
#### propagate

```python
propagate(x, y, cphase, defocus=0, wavelength=1.97e-12, focal_length=1)
```

Calculates the Lorentz image given a specific phase.

Inputs:

x, y: numpy.2darray objects with the x and y coordinates

cphase: the complex phase acquired (complex to allow for attenuation)

defocus: scalar

**Returns**:

  
- `psi_out` - numpy.2darray containing the complex wavefunction in the image plane.

<a name="wsp_tools.lorentzSim.T"></a>
#### T

```python
T(qx, qy, defocus, wavelength)
```

Utility function for propagate(). Microscope transfer function.

<a name="wsp_tools.lorentzSim.chi"></a>
#### chi

```python
chi(qx, qy, defocus, wavelength)
```

Utility function for propagate(). Phase transfer function.

<a name="wsp_tools.lorentzSim.aperture"></a>
#### aperture

```python
aperture(qx, qy)
```

Utility function for propagate(). Circular aperture.

<a name="wsp_tools.lorentzSim.jchessmodel"></a>
#### jchessmodel

```python
jchessmodel(x, y, z=0, **kwargs)
```

Calculates the magnetization of a hopfion based on Jordan Chess' model.

**Returns**:

  
  mx, my, mz - numpy.ndarray objects in the same shape as x and y.

<a name="wsp_tools.beam"></a>
# wsp\_tools.beam

Module to generate spatial modes and beam parameters.

<a name="wsp_tools.beam.E"></a>
#### E

```python
E(T_eV)
```

kinetic energy [eV] -> total energy [J]

<a name="wsp_tools.beam.p"></a>
#### p

```python
p(T_eV)
```

kinetic energy [eV] -> momentum [N s]

<a name="wsp_tools.beam.dB"></a>
#### dB

```python
dB(T_eV)
```

kinetic energy [eV] -> deBroglie wavelength [m]

<a name="wsp_tools.beam.k"></a>
#### k

```python
k(T_eV)
```

kinetic energy [eV] -> wavenumber [m^-1]

<a name="wsp_tools.beam.omega"></a>
#### omega

```python
omega(T_eV)
```

kinetic energy [eV] -> angular frequency [rad s^-1]

<a name="wsp_tools.beam.v_g"></a>
#### v\_g

```python
v_g(T_eV)
```

kinetic energy [eV] -> group velocity [m s^-1]

<a name="wsp_tools.beam.v_p"></a>
#### v\_p

```python
v_p(T_eV)
```

kinetic energy [eV] -> phase velocity [m s^-1]

<a name="wsp_tools.beam.bessel"></a>
#### bessel

```python
bessel(x, y, z=0, l=0, kz0=k(3e5), kperp=0.0005*k(3e5), dkperp=0, N=30)
```

Creates a bessel beam by adding plane waves.

The spectrum is a circle in
k-space (k_z, dkperp + kperp * cos(theta), kperp * sin(theta)).

<a name="wsp_tools.beam.besselPacket"></a>
#### besselPacket

```python
besselPacket(t=0, l=0, kres=2**7, kmin=-3*k(3e5), kmax=3*k(3e5), kz0=k(3e5), kperp=.5*k(3e5), dkperp=0, sig=0.05*k(3e5))
```

Creates a bessel beam by Fourier transforming a Gaussian spectrum.

The spectrum is a gaussian centered on a circle in k-space
(k_z, dkperp + kperp * cos(theta), kperp * sin(theta)).

<a name="wsp_tools.beam.zR"></a>
#### zR

```python
zR(k, w0)
```

Rayleigh range as a function of k, beam waist

<a name="wsp_tools.beam.w"></a>
#### w

```python
w(z, w0, k)
```

Beam waist as a function of z, beam waist at z=0, and k

<a name="wsp_tools.beam.R"></a>
#### R

```python
R(z, w0, k)
```

Radius of curvature as a function of z, w0, k

<a name="wsp_tools.beam.LG"></a>
#### LG

```python
LG(x, y, z=0, l=0, p=0, w_0=2e-6, lam=1.97e-12)
```

Generates a Laguerre-Gauss spatial mode.

**Returns**:

  
  numpy.2darray() = the 2d complex amplitude of a Laguerre-Gauss mode at z.

