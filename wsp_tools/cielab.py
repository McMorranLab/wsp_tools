import numpy as np
from scipy.special import eval_genlaguerre, factorial
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['image.origin'] = 'lower'

def cielab_image(data, alpha):
    rgb_image_components = np.zeros(np.append(data.shape, 4),dtype=np.uint8)
    if alpha == 'intensity':
        value = np.absolute(data)**2 * 255
    else:
        value = np.absolute(data) * 255
    hue = (np.angle(data) + np.pi) / 2
    pi6 = np.pi/6
    def sin2(array, offset):
        return(np.sin(array-offset)**2)
    r = sin2(hue, .15*pi6) + 0.35 * sin2(hue, 3.15 * pi6)
    b = sin2(hue, 4.25 * pi6)
    g = .6*sin2(hue, 2*pi6) + 0.065 * sin2(hue, 5.05 * pi6) + 0.445*b + 0.33*r
    rgb_image_components[:,:,0] = (r * value).astype(np.uint8)
    rgb_image_components[:,:,1] = (g * value).astype(np.uint8)
    rgb_image_components[:,:,2] = (b * value).astype(np.uint8)
    rgb_image_components[:,:,3] = np.full(data.shape, fill_value=255, dtype=np.uint8)
    return(rgb_image_components)

def rgba(mode, cmap = 'uniform', alpha = 'intensity'):
    mode /= np.max(np.abs(mode))
    if cmap == 'uniform':
        out = cielab_image(mode, alpha)
        return(out)
    colormap = plt.cm.ScalarMappable(cmap=cmap)
    out = colormap.to_rgba(np.angle(mode))
    out[:,:,-1] = np.abs(mode)**2/np.max(np.abs(mode)**2)
    return(out)