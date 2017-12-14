### Script to profile two methods, the radial_averaging and the .as_signal

#%matplotlib tk
import pycrystem as pc
import numpy as np
from pycrystem.indexation_generator import *
from hyperspy.components2d import Expression
from pycrystem import ElectronDiffraction

folder = '/home/pc494/Documents/data/Jandrek/'
structure = pc.Structure.from_file(folder+"CsPbBr3_1.cif")
rot_array = np.loadtxt(folder + 'mmm_grid_euler.bin')
rot_list = rot_array.tolist()
rot_list = rot_list[0:500]
edc = pc.ElectronDiffractionCalculator(300, 0.025)
diff_gen = pc.DiffractionLibraryGenerator(edc)
struc_lib = dict()
struc_lib["CsPbBr3"] = (structure, rot_list)
library = diff_gen.get_diffraction_library(struc_lib,
                                            calibration=1.2/128,
                                            reciprocal_radius=1.,
                                            representation='euler')

#scaling is the new problem
radial_sims = {}

_GAUSSIAN2D_EXPR = \
    "intensity * exp(" \
    "-((x-cx)**2 / (2 * sigma ** 2)" \
    " + (y-cy)**2 / (2 * sigma ** 2))" \
")"
@profile
def radial_average(z, center):
    y, x = np.indices(z.shape)
    r = np.sqrt((x - center[1])**2 + (y - center[0])**2)
    r = (r+0.5).astype(np.int)

    tbin = np.bincount(r.ravel(), z.ravel())
    nr = np.bincount(r.ravel())
    averaged = tbin / nr

    return averaged
@profile
def as_signal(z, size, sigma, max_r):
    dp_dat = 0
    l = np.linspace(-max_r, max_r, size)
    x, y = np.meshgrid(l, l)
    coords = z.coordinates[:, :2]
    g = Expression(_GAUSSIAN2D_EXPR, 'Gaussian2D', module='numexpr')
    for (cx, cy), intensity in zip(coords, z.intensities):
        g.intensity.value = intensity
        g.sigma.value = sigma
        g.cx.value = cx
        g.cy.value = cy
        dp_dat += g.function(x, y)

    dp = ElectronDiffraction(dp_dat)
    dp.set_calibration(2*max_r/size)

    return dp

timer = 0
for XZX in library['CsPbBr3']:
    z_temp = as_signal(library['CsPbBr3'][XZX],144,0.1,2.1).data #pull out for profile
    radial_sims[XZX] = radial_average(z_temp,center=[144/2,144/2])
    timer += 1
    if timer%1000==0:
        print(timer)