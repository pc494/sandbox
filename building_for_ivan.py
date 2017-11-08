import pymatgen as pmg
import pyprismatic as pr
import pycrystem as pc
from matplotlib import pyplot as plt
import numpy as np

filepath = '/home/pc494/Documents/code/sandbox/azif.cif'

#struct = pmgcif.CifParser(filepath)
#alpha = struct.get_structures()
try:
    struct = pmg.Structure.from_file(filepath)
except UnicodeDecodeError:
    pass
    
alpha_divider = 3
probe_divider = 3
## Set up our microscope ##
ediff = pc.ElectronDiffractionCalculator(200., 0.025)

## Run our simulations ##
meta_params = {}
meta_params['save4DOutput'] = True
meta_params['save3DOutput'] = False
meta_params['scanWindowXMin'] = 0.495
meta_params['scanWindowXMax'] = 0.50
meta_params['scanWindowYMin'] = 0.495
meta_params['scanWindowYMax'] = 0.50
meta_params['alphaBeamMax'] = 0.024/alpha_divider
meta_params['probeSemiangle'] = 0.02/probe_divider
ediff.calculate_ed_data_dynamic(struct,meta_params)

### Don't mess this up ###
"""
output = pr.fileio.readMRC('PP_output_X1_Y0_FP1.mrc')
output = output.reshape([output.shape[1],output.shape[2]])
output = np.fft.fftshift(output)
plt.figure()
plt.imshow(np.power(output,0.25),cmap='viridis')
plt.draw()
plt.savefig('image.png')
"""

