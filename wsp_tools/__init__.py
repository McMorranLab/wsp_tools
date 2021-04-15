"""
wsp_tools contains utilities for TEM data analysis and presentation.

Features:

* Single Image TIE
* Lorentz simulations - phase calculations, propagation
* spatial mode implementations - LG, Bessel beams, Bessel packets
* basic image processing - high_pass, low_pass, clipping
* a matplotlib.pyplot wrapper
* an implementation of the CIELAB colorspace
* a scipy.constants (CODATA values) wrapper that allows unit scaling (i.e., using nanometers
instead of meters)
"""

import numpy as np
import matplotlib.pyplot as plt
import os, pdoc
from pathlib import Path

from .beam import *
from .cielab import *
from .constants import setUnits
from .image_processing import *
from .lorentzSim import *
from .pyplotwrapper import *
from .sitie import *

def docs(outdir = "."):
	"""Auto-generate documentation for wsp-tools in html.

	**Parameters**

	* **outdir** : _string_ <br />
	The directory to write the output documentation. <br />
	Default is "./".
	"""
	modules = ['wsp_tools']
	context = pdoc.Context()
	modules = [pdoc.Module(mod, context=context)
			   for mod in modules]
	pdoc.link_inheritance(context)

	def recursive_htmls(mod):
		yield mod.url(), mod.html()
		for submod in mod.submodules():
			yield from recursive_htmls(submod)

	for mod in modules:
		for module_url, html in recursive_htmls(mod):
			output_url = Path(outdir).expanduser().joinpath(module_url)
			if not Path(output_url).parent.exists():
				Path(output_url).parent.mkdir(parents=True)
			with open(output_url, "w+") as f:
				f.write(html)
	print("Documentation for wsp-tools written to: \n{}:".format(Path(outdir).joinpath(modules[0].url())))
