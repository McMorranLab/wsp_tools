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
from deprecated import deprecated

from .beam import *
from .cielab import *
from .constants import setUnits
from .image_processing import *
from .pyplotwrapper import *
from .sitie import *
from .ltem import *

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
