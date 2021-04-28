# wsp_tools
wsp_tools contains utilities for TEM data analysis and presentation.

Features:

* Single Image TIE
* Lorentz simulations
* spatial modes
* a `matplotlib.pyplot` wrapper
* an implementation of the CIELAB colorspace
* a `scipy.constants` wrapper that allows arbitrary units (i.e., using nanometers instead of meters)

## Installation
```Bash
pip install -e git+https://github.com/McMorranLab/wsp_tools#egg=wsp-tools
```

Installing in editable mode `-e` ensures that pip records the install url, with the correct commit.

## Documentation

As of version 1.0.94, wsp_tools includes a helper function `docs(outdir=".")`, which uses ```pdoc3``` to auto-generate html documentation and save it to an outdir of your choice, thus avoiding the need to clone this repository just to see the documentation. For example:

```Python
import wsp_tools
wsp_tools.docs("~/Desktop") # generates docs in the folder ~/Desktop/wsp_tools
```

Documentation is auto-generated by ```pdoc``` and is included in two forms:

1. [A single markdown file](docs/docsmd/docs.md) (for readability on GitHub)
2. [As html](docs/wsp_tools/index.html) (easier to navigate, but you have to clone the repository to render them).

```Bash
git clone https://github.com/McMorranLab/wsp_tools
open wsp_tools/docs/wsp_tools/index.html
```

## Tests

Tests are split into two subdirectories:

1. `tests`
	These are typical unit tests, that assert that functions return the right shape, beam parameters return the right values, etc.
2. `devtests`
	These are tests of the actual functionality, that require a trained eye to evaluate. For example, a function `test_bessel()` will generate a bessel beam using `wsp_tools.bessel()`, but rather than asserting a unit test, will just plot the bessel beam so the developer can manually assert whether it looks correct.

The rationale for `devtests` is that this package is math-heavy, so it's highly possible for the code to run fine, but be wrong. The easiest way to test for this is to check base cases where the developer knows what to look for.

Tests can be performed by

```Bash
cd tests
pytest
```

and

```Bash
cd devtests
bash devtests.sh
```

or
```Bash
python test_beam.py
python test_cielab.py
... etc
```
