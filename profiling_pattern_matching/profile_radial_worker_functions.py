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

def radial_average(z, center):
    y, x = np.indices(z.shape)
    r = np.sqrt((x - center[1])**2 + (y - center[0])**2)
    r = (r+0.5).astype(np.int)

    tbin = np.bincount(r.ravel(), z.ravel())
    nr = np.bincount(r.ravel())
    averaged = tbin / nr

    return averaged


def which_edge(l,d,d_up):
    if (d - l[d_up-1]) > (l[d_up]-d):
        return d_up
    else:
        return d_up-1
def as_pure_peaks(z,size,max_r):
    l = np.linspace(-max_r, max_r, size)
    coords = z.coordinates[:, :2]
    signal = np.zeros([size,size])
    for i in np.arange(coords.shape[0]):
        x,y = coords[i,0],coords[i,1]
        x_up,y_up = np.sum(l < x),np.sum(l < y) # when x > l we have overshot slightly, 
        x_num,y_num = which_edge(l,x,x_up),which_edge(l,y,y_up)
        ## next fix the intensity
        signal[x_num,y_num] += z.intensities[i]
    dp = ElectronDiffraction(signal)
    dp.set_calibration(2*max_r/size)
    
    return dp

    

timer = 0
for XZX in library['CsPbBr3']:
    z_temp = as_pure_peaks(library['CsPbBr3'][XZX],144,2.1).data #pull out for profile
    radial_sims[XZX] = radial_average(z_temp,center=[144/2,144/2])
    timer += 1
    if timer%1000==0:
        print(timer)