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

"""Wrapper for the scipy.constants module that allows unit scaling.

To see all available constants, `print(wsp_tools.constants.__all__)`.

To set units:

```python
import wsp_tools as wt
wt.setUnits(meter = 1e-3) # set km as base unit for length

print(wt.constants.c) # outputs 299792.458
```

Note that all other modules should update automatically as well.
"""

import scipy.constants as sc

__all__ = ['s','m','kg','A','K','mol','cd','F','J','C','W','eV','c',
			'mu_0','eps_0','h','hbar','h_eV','hbar_eV','G','e','R','alpha',
			'N_A','k_B','k_B_eV','sigma','Wien','Rydberg','m_e','m_p',
			'm_n','pi','setUnits']

def setUnits(second=1, meter=1, kilogram=1, Amp=1, Kelvin=1, mole=1, candela=1):
	"""Sets the units across the wsp_tools module.

	i.e. setUnits(meter = 1000) sets the millimeter as the base unit for length.

	**Parameters**

	* **second** : _number, optional_ <br />
	The SI base unit for time. <br />
	Default is `second = 1`.

	* **meter** : _number, optional_ <br />
	The SI base unit for length. <br />
	Default is `meter = 1`.

	* **kilogram** : _number, optional_ <br />
	The SI base unit for mass. <br />
	Default is `kilogram = 1`.

	* **Amp** : _number, optional_ <br />
	The SI base unit for current. <br />
	Default is `Amp = 1`.

	* **Kelvin** : _number, optional_ <br />
	The SI base unit for temperature. <br />
	Default is `Kelvin = 1`.

	* **mole** : _number, optional_ <br />
	The SI base unit for amount of substance. <br />
	Default is `mole = 1`.

	* **candela** : _number, optional_ <br />
	The SI base unit for luminous intensity. <br />
	Default is `candela = 1`.
	"""
	global s,m,kg,A,K,mol,cd
	s,m,kg,A,K,mol,cd = second, meter, kilogram, Amp, Kelvin, mole, candela
	setConsts(s, m, kg, A, K, mol, cd)

def setConsts(s, m, kg, A, K, mol, cd):
	"""Utility function for `setUnits()`."""
	global F,J,C,W,eV,c,mu_0,eps_0,h,hbar,h_eV,hbar_eV,G,e,R,alpha,N_A,k_B
	global k_B_eV,sigma,Wien,Rydberg,m_e,m_p,m_n,pi
	F = s**4 * A**2 / m**2 / kg
	J = kg * m**2 / s**2
	C = A * s
	W = kg * m**2 / s**3
	eV = sc.e * J

	c       = sc.c          * m / s
	mu_0    = sc.mu_0       * m * kg / s**2 /A**2
	eps_0   = sc.epsilon_0  * F / m
	h       = sc.h          * J * s
	hbar    = sc.hbar       * J * s
	h_eV    = sc.h          * s * J / eV
	hbar_eV = sc.hbar       * s * J / eV
	G       = sc.G          * m**3 / kg / s**2
	g       = sc.g          * m / s**2
	e       = sc.e          * C
	R       = sc.R          * kg / s**2 / K / mol
	alpha   = sc.alpha
	N_A     = sc.N_A        / mol
	k_B     = sc.k          * J / K
	k_B_eV  = sc.k          * J / eV / K
	sigma   = sc.sigma      * W / m**2 / K**4
	Wien    = sc.Wien       * m * K
	Rydberg = sc.Rydberg    / m
	m_e     = sc.m_e        * kg
	m_p     = sc.m_p        * kg
	m_n     = sc.m_n        * kg
	pi      = sc.pi

setUnits()
