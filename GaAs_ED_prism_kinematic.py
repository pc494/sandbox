import pycrystem as pc
import pymatgen as pmg
from matplotlib import pyplot as plt
import numpy as np
from pycrystem.utils.pyprismatic_io_utils import *
from matplotlib import pyplot as plt
import pyprismatic as pr

## Set up our structure ##
Ga = pmg.Element("Ga")
As = pmg.Element("As")
lattice = pmg.Lattice.cubic(5.65)
structure = pmg.Structure.from_spacegroup("F-43m",lattice, [Ga,As], [[0, 0, 0],[0.5,0.5,0.5]])
structure.make_supercell([4,4,1]) #square to remain compatible with kinematic

## Set up our microscope ##
ediff = pc.ElectronDiffractionCalculator(200., 0.025)

## Run our simulations ##
meta_params = {}
meta_params['save4DOutput'] = True
meta_params['save3DOutput'] = False
meta_params['scanWindowXMin'] = 0.49
meta_params['scanWindowXMax'] = 0.51
meta_params['scanWindowYMin'] = 0.49
meta_params['scanWindowYMax'] = 0.51
#meta_params['realspacePixelSizeX']=0.1
#meta_params['realspacePixelSizeY']=0.1

x_tiles = 5
y_tiles = 0

#ediff.calculate_ed_data_dynamic(structure,meta_params)

def create_horizontal(y_line,x_len):
    image = np.squeeze(pr.fileio.readMRC('PP_output_X0_Y'+str(y_line)+'_FP1.mrc'))
    for x_cord in np.arange(1,x_len):
        file_string = 'PP_output_X' + str(x_cord) + '_Y' + str(y_line) + '_FP1.mrc'
        image = np.concatenate([image,np.squeeze(pr.fileio.readMRC(file_string))],axis=1)
    return image

image = create_horizontal(0,4)
print(image.shape)
for y in [1,2,3]:
    image = np.concatenate([image,create_horizontal(y,4)],axis=0)
    print(image.shape)

plt.imshow(image)
plt.show()
"""
output = pr.fileio.readMRC('PP_output_X1_Y0_FP1.mrc')
output3 = pr.fileio.readMRC('PP_output_X1_Y0_FP1.mrc')
output4 = pr.fileio.readMRC('PP_output_X1_Y1_FP1.mrc')
output = output.reshape([output.shape[1],output.shape[2]])
output2 = output.reshape([output2.shape[1],output2.shape[2]])
output = np.concatenate([output,output2],axis=1)
plt.imshow(output)
plt.show()
"""
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

