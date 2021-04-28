import wsp_tools as wt
from wsp_tools import plt, np
import pickle

# %%
def load_data():
	with open("001_parallel.pkl", "rb") as f:
		return(pickle.load(f))

file = load_data()
file['data'] = np.array(file['data'])
file['coords'] = np.array(file['coords'])

img = wt.lorentz(file)

# %%
img.data.clip_data().high_pass().shift_pos()
img.sitie(defocus=1e-3)

def test_phase():
	assert(len(img.phase.shape) == 2)

def test_B():
	assert(len(img.Bx.shape) == 2)

def test_units():
	assert(img.xUnit == 'm')
	assert(np.isclose(img.dx,2.443537348881364e-09))
