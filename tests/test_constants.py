import wsp_tools as wt
from wsp_tools import constants as _

# %%
def test_setUnits():
	assert(wt.np.isclose(_.c,299792458.0))
	wt.setUnits(meter=2)
	assert(wt.np.isclose(_.c,599584916.0))
