# wsp_tools
Utilities for TEM data analysis and simulation.

Features:

* Single Image TIE
* Lorentz TEM simulations
* spatial modes
* a `matplotlib.pyplot` wrapper
* an implementation of the CIELAB colorspace
* a `scipy.constants` wrapper that allows arbitrary units (i.e., using nanometers instead of meters)

## Installation

### Pip
As of version 1.3.1, `wsp_tools` is available on PyPi, and can be installed with `pip`:

```Bash
python -m pip install wsp-tools
```

### GitHub
`wsp_tools` can also be installed directly from GitHub via
```Bash
pip install -e git+https://github.com/McMorranLab/wsp_tools#egg=wsp-tools
```

Installing in editable mode `-e` ensures that `pip freeze` records the install url, not just the version.

You can specify the commit or tag to install, via

```Bash
pip install -e git+https://github.com/McMorranLab/wsp_tools@{commit id or tag number}#egg=wsp_tools
```

## Documentation

As of version 1.3.1, documentation is available at [https://mcmorranlab.github.io/wsp_tools/](https://mcmorranlab.github.io/wsp_tools/).

As of version 1.0.94, wsp_tools includes a helper function `docs(outdir=".")`, which uses ```pdoc3``` to auto-generate html documentation and save it to an outdir of your choice, thus avoiding the need to clone this repository just to see the documentation. For example:

```Python
import wsp_tools
wsp_tools.docs("~/Desktop") # generates docs in the folder ~/Desktop/wsp_tools
```

## Tests

Tests are split into two subdirectories:

1. `tests`
	These are typical unit tests, that assert that functions return the right shape, beam parameters return the right values, etc. Run with `pytest`.
2. `devtests`
	These are tests of the actual functionality, that require a trained eye to evaluate. For example, a function `test_bessel()` will generate a bessel beam using `wsp_tools.bessel()`, but rather than asserting a unit test, will just plot the bessel beam so the developer can manually assert whether it looks correct. Run as normal `.py` scripts.

The rationale for `devtests` is that this package is math-heavy, so it's highly possible for the code to run fine, but be wrong. The easiest way to test for this is to check base cases where the developer knows what to look for.
