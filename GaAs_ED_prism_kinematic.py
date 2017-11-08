#!/usr/bin/env python3

""" Usage:
     GaAs_ED_prism_kinematic.py <alpha_divider> <probe_divider> <xy_size> <z_size>
"""

from docopt import docopt
import pycrystem as pc
import pymatgen as pmg
from matplotlib import pyplot as plt
import numpy as np
from pycrystem.utils.pyprismatic_io_utils import *
from matplotlib import pyplot as plt
import pyprismatic as pr

args = docopt(__doc__)

alpha_divider = float(args["<alpha_divider>"]) #3 or nearer to 1 than that
probe_divider = float(args["<probe_divider>"]) #3 or nearer to 1 than that
xy_size = float(args["<xy_size>"])
z_size  = float(args["<z_size>"])


## Set up our structure ##
Ga = pmg.Element("Ga")
As = pmg.Element("As")
lattice = pmg.Lattice.cubic(5.65)
structure = pmg.Structure.from_spacegroup("F-43m",lattice, [Ga,As], [[0, 0, 0],[0.5,0.5,0.5]])
structure.make_supercell([xy_size,xy_size,z_size]) #square to remain compatible with kinematic

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
ediff.calculate_ed_data_dynamic(structure,meta_params)

output = pr.fileio.readMRC('PP_output_X1_Y0_FP1.mrc')
output = output.reshape([output.shape[1],output.shape[2]])
output = np.fft.fftshift(output)
plt.figure()
plt.imshow(np.power(output,0.25),cmap='viridis')
plt.draw()
plt.savefig('image.png')
"""
# The line above has been excuted, the output is stored in an .mrc file
structure = pmg.Structure.from_spacegroup("F-43m",lattice, [Ga,As], [[0, 0, 0],[0.5,0.5,0.5]]) 
structure.make_supercell([4,4,1]) #no depth needed for kinematic
dd_kinematic = ediff.calculate_ed_data_kinematic(structure,reciprocal_radius=2.5)

## Post Processing ##

ddd_buffer = import_pyprismatic_data(meta)
pltable_dd_d = np.squeeze(np.sum(ddd_buffer,axis=2))
plt.imshow(pltable_dd_d)
plt.draw()

pltable_dd_k = pc.ElectronDiffraction(dd_kinematic.as_signal(512,0.05,2.5).data)
pltable_dd_k.plot()
plt.draw()
"""


#plt.show() #to hold figures at the end of the script
