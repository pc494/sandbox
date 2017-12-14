### Avoiding a notebook as it won't play nice with tqdm

#%matplotlib tk
import pycrystem as pc
import numpy as np
import hyperspy.api as hs
from pycrystem.indexation_generator import *
from pycrystem.utils.expt_utils import *
from matplotlib import pyplot as plt
from pycrystem import ElectronDiffraction

folder = '/home/pc494/Documents/data/Jandrek/'

dp = pc.load(folder+'QD-03.blo')

dp.apply_affine_transformation(D = np.array([[0.99, 0.00, 0.00],
                                             [0.00, 0.69, 0.00],
                                             [0.00, 0.00, 1.00]]),show_progressbar=False)
dp.set_calibration(0.032)

## Take radial averages - note the non-square error might be okay
radial_averaged_data = dp.map(radial_average,center=np.asarray([144/2,144/2]),cython=False,inplace=False, show_progressbar=False)
print(radial_averaged_data)

dp.plot()
radial_averaged_data_signal = hs.signals.Signal1D(radial_averaged_data)
radial_averaged_data_signal.plot()

structure = pc.Structure.from_file(folder+"CsPbBr3_1.cif")
rot_array = np.loadtxt(folder + 'mmm_grid_euler.bin')
rot_list = rot_array.tolist()
edc = pc.ElectronDiffractionCalculator(300, 0.025)
diff_gen = pc.DiffractionLibraryGenerator(edc)
struc_lib = dict()
struc_lib["CsPbBr3"] = (structure, rot_list)
library = diff_gen.get_diffraction_library(struc_lib,
                                            calibration=1.2/128,
                                            reciprocal_radius=1.,
                                            representation='euler')

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


pattern = library['CsPbBr3'][(1.236036,1.003564,3.563619)]
plt.show()
#scaling is the new problem
radial_sims = {}
#(radial_average(pattern.as_signal(144,0.1,2.1).data,center=[144/2,144/2],cython=False))
timer = 0
for XZX in library['CsPbBr3']:
    peaks_2D_data = as_pure_peaks(library['CsPbBr3'][XZX],144,2.1).data
    radial_sims[XZX] = radial_average(peaks_2D_data,center=[144/2,144/2],cython=False)
    timer += 1
    if timer%1000==0:
        print(timer)