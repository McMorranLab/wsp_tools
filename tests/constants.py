import os, sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import wsp_tools as wt
from wsp_tools import constants as _

# %%
print(_.c)
wt.setUnits(meter=2)
print(_.c, "meter set to 2")
