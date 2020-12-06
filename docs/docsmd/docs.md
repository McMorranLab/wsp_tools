---
description: |
    API documentation for modules: wsp_tools, wsp_tools.beam, wsp_tools.cielab, wsp_tools.constants, wsp_tools.lorentzSim, wsp_tools.pyplotwrapper, wsp_tools.sitie.

lang: en

classoption: oneside
geometry: margin=1in
papersize: a4

linkcolor: blue
links-as-notes: true
...


    
# Module `wsp_tools` {#wsp_tools}

wsp_tools contains utilities for TEM data analysis and presentation.

Features:

* Single Image TIE
* Lorentz simulations
* spatial mode implementations
* a matplotlib.pyplot wrapper
* an implementation of the cielab colorspace
* a scipy.constants wrapper that allows unit scaling (i.e., using nanometers
instead of meters)


    
## Sub-modules

* [wsp_tools.beam](#wsp_tools.beam)
* [wsp_tools.cielab](#wsp_tools.cielab)
* [wsp_tools.constants](#wsp_tools.constants)
* [wsp_tools.lorentzSim](#wsp_tools.lorentzSim)
* [wsp_tools.pyplotwrapper](#wsp_tools.pyplotwrapper)
* [wsp_tools.sitie](#wsp_tools.sitie)






    
# Module `wsp_tools.beam` {#wsp_tools.beam}

Module to generate spatial modes and beam parameters.




    
## Functions


    
### Function `E` {#wsp_tools.beam.E}




>     def E(
>         T_eV
>     )


kinetic energy [eV] -> total energy [J]

    
### Function `LG` {#wsp_tools.beam.LG}




>     def LG(
>         x,
>         y,
>         z=0,
>         l=0,
>         p=0,
>         w_0=2e-06,
>         lam=1.97e-12
>     )


Generates a Laguerre-Gauss spatial mode.

Returns:

numpy.2darray() = the 2d complex amplitude of a Laguerre-Gauss mode at z.

    
### Function `R` {#wsp_tools.beam.R}




>     def R(
>         z,
>         w0,
>         k
>     )


Radius of curvature as a function of z, w0, k

    
### Function `bessel` {#wsp_tools.beam.bessel}




>     def bessel(
>         x,
>         y,
>         z=0,
>         l=0,
>         kz0=3191460985702.0464,
>         kperp=1595730492.8510232,
>         dkperp=0,
>         N=30
>     )


Creates a bessel beam by adding plane waves.

The spectrum is a circle in
k-space (k_z, dkperp + kperp * cos(theta), kperp * sin(theta)).

    
### Function `besselPacket` {#wsp_tools.beam.besselPacket}




>     def besselPacket(
>         t=0,
>         l=0,
>         kres=128,
>         kmin=-9574382957106.139,
>         kmax=9574382957106.139,
>         kz0=3191460985702.0464,
>         kperp=1595730492851.0232,
>         dkperp=0,
>         sig=159573049285.10233
>     )


Creates a bessel beam by Fourier transforming a Gaussian spectrum.

The spectrum is a gaussian centered on a circle in k-space
(k_z, dkperp + kperp * cos(theta), kperp * sin(theta)).

    
### Function `dB` {#wsp_tools.beam.dB}




>     def dB(
>         T_eV
>     )


kinetic energy [eV] -> deBroglie wavelength [m]

    
### Function `k` {#wsp_tools.beam.k}




>     def k(
>         T_eV
>     )


kinetic energy [eV] -> wavenumber [m^-1]

    
### Function `omega` {#wsp_tools.beam.omega}




>     def omega(
>         T_eV
>     )


kinetic energy [eV] -> angular frequency [rad s^-1]

    
### Function `p` {#wsp_tools.beam.p}




>     def p(
>         T_eV
>     )


kinetic energy [eV] -> momentum [N s]

    
### Function `v_g` {#wsp_tools.beam.v_g}




>     def v_g(
>         T_eV
>     )


kinetic energy [eV] -> group velocity [m s^-1]

    
### Function `v_p` {#wsp_tools.beam.v_p}




>     def v_p(
>         T_eV
>     )


kinetic energy [eV] -> phase velocity [m s^-1]

    
### Function `w` {#wsp_tools.beam.w}




>     def w(
>         z,
>         w0,
>         k
>     )


Beam waist as a function of z, beam waist at z=0, and k

    
### Function `zR` {#wsp_tools.beam.zR}




>     def zR(
>         k,
>         w0
>     )


Rayleigh range as a function of k, beam waist




    
# Module `wsp_tools.cielab` {#wsp_tools.cielab}

Module to generate rgba data from numpy.2darrays.




    
## Functions


    
### Function `cielab_image` {#wsp_tools.cielab.cielab_image}




>     def cielab_image(
>         data,
>         alpha='intensity'
>     )


Converts a 2d complex array to rgba data based on the cielab color space.

Input:

numpy.2darray(). dtype can be real or complex - if real, the output color
will be constant. If complex, the color will reflect the complex angle.

Returns:

numpy.ndarray() with shape [data.shape[0], data.shape[1], 4]

    
### Function `rgba` {#wsp_tools.cielab.rgba}




>     def rgba(
>         mode,
>         cmap='uniform',
>         alpha='intensity'
>     )


Converts a 2d complex array to rgba data.

Input:

dtype can be real or complex - if real, the output color
will be constant. If complex, the color will reflect the complex angle.

If cmap = 'uniform', the cielab color space will be used.

Alpha can be either intensity or amplitude.

Output:

numpy.ndarray() with shape [data.shape[0], data.shape[1], 4].

    
### Function `whatIsC` {#wsp_tools.cielab.whatIsC}




>     def whatIsC()


Used for testing the constants functionality of the module.




    
# Module `wsp_tools.constants` {#wsp_tools.constants}

Wrapper for the scipy.constants module that allows unit scaling.

To see all available constants, print(wsp_tools.constants.__all__).

To set units:

import wsp_tools as wt
wt.setUnits(meter = 1e-3) # set km as base unit for length
from wsp_tools.constants import *

print(c) # outputs 299792.458

Note that all other modules should update automatically as well.







    
# Module `wsp_tools.lorentzSim` {#wsp_tools.lorentzSim}






    
## Functions


    
### Function `T` {#wsp_tools.lorentzSim.T}




>     def T(
>         qx,
>         qy,
>         defocus,
>         wavelength
>     )


Utility function for propagate(). Microscope transfer function.

    
### Function `abphase2d` {#wsp_tools.lorentzSim.abphase2d}




>     def abphase2d(
>         mx,
>         my,
>         mz,
>         Lx=1e-06,
>         Ly=1e-06,
>         p=array([0, 0, 1]),
>         t=6e-08
>     )


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

Returns:

phi = numpy.2darray() the phase acquired.

    
### Function `aperture` {#wsp_tools.lorentzSim.aperture}




>     def aperture(
>         qx,
>         qy
>     )


Utility function for propagate(). Circular aperture.

    
### Function `chi` {#wsp_tools.lorentzSim.chi}




>     def chi(
>         qx,
>         qy,
>         defocus,
>         wavelength
>     )


Utility function for propagate(). Phase transfer function.

    
### Function `jchessmodel` {#wsp_tools.lorentzSim.jchessmodel}




>     def jchessmodel(
>         x,
>         y,
>         z=0,
>         **kwargs
>     )


Calculates the magnetization of a hopfion based on Jordan Chess' model.

Returns:

mx, my, mz - numpy.ndarray objects in the same shape as x and y.

    
### Function `propagate` {#wsp_tools.lorentzSim.propagate}




>     def propagate(
>         x,
>         y,
>         cphase,
>         defocus=0,
>         wavelength=1.97e-12,
>         focal_length=1
>     )


Calculates the Lorentz image given a specific phase.

Inputs:

x, y: numpy.2darray objects with the x and y coordinates

cphase: the complex phase acquired (complex to allow for attenuation)

defocus: scalar

Returns:

psi_out: numpy.2darray containing the complex wavefunction in the image plane.




    
# Module `wsp_tools.pyplotwrapper` {#wsp_tools.pyplotwrapper}

A wrapper for matplotlib.pyplot containing common plotting routines.




    
## Functions


    
### Function `subplots` {#wsp_tools.pyplotwrapper.subplots}




>     def subplots(
>         rc=11,
>         sharex=False,
>         sharey=False,
>         squeeze=False,
>         subplot_kw=None,
>         gridspec_kw=None,
>         **fig_kw
>     )


Creates a fig, ax instance but replaces ax with singleAx.

Behaves identically to matplotlib.pyplot.subplots(), but replaces each
matplotlib.axes.Axes object with a wsp_tools.pyplotwrapper.singleAx
object.

Each wsp_tools.pyplotwrapper.singleAx object in turn behaves just like a
normal Axes object, but with added methods.


    
## Classes


    
### Class `singleAx` {#wsp_tools.pyplotwrapper.singleAx}




>     class singleAx(
>         ax,
>         title='',
>         xlabel='',
>         ylabel=''
>     )


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

More commonly, this class is returned by ```wsp_tools.pyplotwrapper.subplots```.







    
#### Methods


    
##### Method `imshow` {#wsp_tools.pyplotwrapper.singleAx.imshow}




>     def imshow(
>         self,
>         data,
>         step=1,
>         **kwargs
>     )


Imshows the data.

If a window has been applied to the plot before this is called, it
will be used.

data should be a numpy.2darray object.
step: data[::step,::step] will be shown.

    
##### Method `inset` {#wsp_tools.pyplotwrapper.singleAx.inset}




>     def inset(
>         self,
>         window=None,
>         color='white',
>         **kwargs
>     )


Plots a square box with vertices defined by window.

Window = (xmin, xmax, ymin, ymax).

Default color is white. Takes all the same kwargs as
matplotlib.pyplot.plot().

    
##### Method `prePlot` {#wsp_tools.pyplotwrapper.singleAx.prePlot}




>     def prePlot(
>         self,
>         data,
>         step=1
>     )


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

Returns:

1. xout - numpy.1darray: windowed x-axis
2. yout - numpy.1darray: windowed y-axis
3. dout - numpy.2darray: windowed data

Additionally, it sets the ax element's extent.

    
##### Method `quiver` {#wsp_tools.pyplotwrapper.singleAx.quiver}




>     def quiver(
>         self,
>         data,
>         step=1,
>         **kwargs
>     )


Shows a quiver plot of complex data.

data should be a complex numpy.2darray object.
        If real, the y component will just be zero everywhere.
step: data[::step,::step] will be shown.

    
##### Method `rgba` {#wsp_tools.pyplotwrapper.singleAx.rgba}




>     def rgba(
>         self,
>         data,
>         step=1,
>         alpha='intensity',
>         cmap='uniform',
>         **kwargs
>     )


Shows an rgba interpretation of complex data.

Color represents complex angle, and brightness represents either
amplitude or intensity.

data should be a complex numpy.2darray object.
        If real, the color will be uniform and only brightness will vary.

alpha can represent either 'intensity' (default) or 'amplitude'.

    
##### Method `setAxes` {#wsp_tools.pyplotwrapper.singleAx.setAxes}




>     def setAxes(
>         self,
>         x,
>         y,
>         window=None
>     )


Sets the x and y axes of the singleAx object, and can apply a window.

x, y should be numpy.1darray objects.

Window should be wrt the values of x and y, not their
indices.
If window is omitted, the full data and x-y axes will be used.

    
##### Method `setWindow` {#wsp_tools.pyplotwrapper.singleAx.setWindow}




>     def setWindow(
>         self,
>         window=(0, 100, 0, 100)
>     )


Applies a window to the singleAx object.

Can be applied before or after setAxes - whichever was applied last
will be used.

window = (xmin, xmax, ymin, ymax) wrt the values of the axes, not
their indices.

    
##### Method `set_title` {#wsp_tools.pyplotwrapper.singleAx.set_title}




>     def set_title(
>         self,
>         title='',
>         **kwargs
>     )


Sets the title of the plot.

Takes all the same args as matplotlib.axes.Axes.set_title

    
##### Method `set_xlabel` {#wsp_tools.pyplotwrapper.singleAx.set_xlabel}




>     def set_xlabel(
>         self,
>         xlabel='',
>         **kwargs
>     )


Sets the xlabel of the plot.

Takes all the same args as matplotlib.axes.Axes.set_xlabel

    
##### Method `set_xytitle` {#wsp_tools.pyplotwrapper.singleAx.set_xytitle}




>     def set_xytitle(
>         self,
>         xlabel='',
>         ylabel='',
>         title='',
>         **kwargs
>     )


Set the xlabel, ylabel, and title at the same time.

Sets all three even if not all are given.
Takes all the same kwargs as matplotlib.axes.Axes.set_xlabel,
matplotlib.axes.Axes.set_ylabel, matplotlib.axes.Axes.set_title.
Whatever you input will be applied to all three.

For individual control, use singleAx.set_xlabel, singleAx.set_ylabel,
or singleAx.set_title.

    
##### Method `set_ylabel` {#wsp_tools.pyplotwrapper.singleAx.set_ylabel}




>     def set_ylabel(
>         self,
>         ylabel='',
>         **kwargs
>     )


Sets the ylabel of the plot.

Takes all the same kwargs as matplotlib.axes.Axes.set_ylabel



    
# Module `wsp_tools.sitie` {#wsp_tools.sitie}

Contains utilities for reconstructing phase and magnetization from Lorentz images.

The most common use case is to generate a lorentz object from a ```.dm3``` file.
Then one can analyze using high_pass(), sitie(), crop_pixel_counts(), etc.

Example:

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




    
## Functions


    
### Function `B_from_phase` {#wsp_tools.sitie.B_from_phase}




>     def B_from_phase(
>         phase,
>         thickness=1
>     )


Calculates the transverse B-field that would impart a specific phase.

    
### Function `SITIE` {#wsp_tools.sitie.SITIE}




>     def SITIE(
>         image,
>         defocus,
>         pixel_size,
>         wavelength=1.97e-12
>     )


Reconstruct the phase from a defocussed image.

    
### Function `blur` {#wsp_tools.sitie.blur}




>     def blur(
>         image,
>         sigma=5,
>         mode='wrap',
>         cval=0.0
>     )


Applies a Gaussian filter to the image.

    
### Function `crop_pixel_counts` {#wsp_tools.sitie.crop_pixel_counts}




>     def crop_pixel_counts(
>         image,
>         sigma=10
>     )


Crops the pixel counts to avg +/- sigma*std.

    
### Function `high_pass` {#wsp_tools.sitie.high_pass}




>     def high_pass(
>         image,
>         sigma=50
>     )


Applies a high-pass filter to the image data.

    
### Function `inverse_laplacian` {#wsp_tools.sitie.inverse_laplacian}




>     def inverse_laplacian(
>         f,
>         pixel_size
>     )




    
### Function `low_pass` {#wsp_tools.sitie.low_pass}




>     def low_pass(
>         image,
>         sigma=50
>     )


Applies a low-pass filter to the image data.

    
### Function `sitie_RHS` {#wsp_tools.sitie.sitie_RHS}




>     def sitie_RHS(
>         I,
>         defocus,
>         wavelength=1.9687489006848795e-12
>     )





    
## Classes


    
### Class `lorentz` {#wsp_tools.sitie.lorentz}




>     class lorentz(
>         dm3file
>     )


Class that contains sitie information about a lorentz image.

Input:

* dm3file: a dictionary-like object with the following keys:
        * data: numpy.2darray() containing the electron counts
        * pixelSize: [scalar, scalar] containing the x and y pixel sizes
        * pixelUnit: [string, string] containing the unit of the pixel sizes







    
#### Methods


    
##### Method `blur` {#wsp_tools.sitie.lorentz.blur}




>     def blur(
>         self,
>         sigma=5,
>         mode='wrap',
>         cval=0.0
>     )


Applies a Gaussian blur to the image data.

    
##### Method `crop_pixel_counts` {#wsp_tools.sitie.lorentz.crop_pixel_counts}




>     def crop_pixel_counts(
>         self,
>         sigma=10
>     )


Crops any pixel counts that are higher or lower than some std from avg.

Sets those pixels to avg +/- sigma*std.

    
##### Method `high_pass` {#wsp_tools.sitie.lorentz.high_pass}




>     def high_pass(
>         self,
>         sigma=20
>     )


Applies a high-pass filter to the image data.

    
##### Method `low_pass` {#wsp_tools.sitie.lorentz.low_pass}




>     def low_pass(
>         self,
>         sigma=50
>     )


Applies a low-pass filter to the image data.

    
##### Method `preview` {#wsp_tools.sitie.lorentz.preview}




>     def preview(
>         self,
>         window=(0, -1, 0, -1)
>     )


Preview the image.

window = (xmin, xmax, ymin, ymax)

    
##### Method `reset` {#wsp_tools.sitie.lorentz.reset}




>     def reset(
>         self
>     )


Resets data to the rawData.

    
##### Method `saveMeta` {#wsp_tools.sitie.lorentz.saveMeta}




>     def saveMeta(
>         self,
>         outdir='',
>         outname='metadata.json'
>     )


Save the metadata of the lorentz object to a file.

    
##### Method `sitie` {#wsp_tools.sitie.lorentz.sitie}




>     def sitie(
>         self,
>         defocus,
>         wavelength=1.96e-12
>     )


Carries out phase and B-field reconstruction.

Assigns phase, Bx, and By attributes.


-----
Generated by *pdoc* 0.9.2 (<https://pdoc3.github.io>).
