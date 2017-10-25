#first commit

from ase import atoms  # may also be worth checking for binneed to check for the binaries perhaps
from pyqstem import PyQSTEM
import pymatgen as pmg
import numpy as np
from pymatgen.io import ase as pymatgenase
from matplotlib import pyplot as plt
import hyperspy.api as hs
from pyqstem.imaging import CTF
from pyqstem.wave import view


si = pmg.Element("Si")
lattice = pmg.Lattice.cubic(5.431)
structure = pmg.Structure.from_spacegroup("Fd-3m",lattice, [si], [[0, 0, 0]])
structure.make_supercell([20,20,10])

accelerating_voltage = 300

qstem = PyQSTEM('TEM') #initialise a TEM object
qstem.set_atoms(pymatgenase.AseAtomsAdaptor.get_atoms(structure)) #this does a pymatgen ---> conversion
qstem.build_wave('plane',accelerating_voltage,(144,144)) ## where do these hardwired numbers come from
num_slices=10 ### what should this number be?
qstem.build_potential(num_slices)
qstem.run() ## this now means we have a wave that has passed through our sample
wave=qstem.get_wave()
ctf = CTF(defocus=50,Cs=5*10**4,focal_spread=30)
#wave=wave.apply_ctf(ctf)

#view(wave,'diffraction pattern')
########################################

def convert_wave_to_diffraction_pattern_array(wave):
    array=wave.array
    img = (np.abs(np.fft.fftshift(np.fft.fft2(array))))**2
    #img = np.log(np.abs(np.fft.fftshift(np.fft.fft2(array))))**2
    extent=wave.get_reciprocal_extent()
    img=img[None,:,:]
    img=img[-1,:,:].T
    img[72][72] = 0 
    return img,extent 

img,img_extent = convert_wave_to_diffraction_pattern_array(wave)
fig, ax = plt.subplots(figsize=(6,6))
imshow = ax.imshow(img)
plt.figure()

#high_values_flags = img < 0.5  # Where values are low
#img[high_values_flags] = 0  # All low values set to 0
spaced_values = np.linspace(img_extent[0],img_extent[1],num=img.shape[0])
intersection_mesh = np.asarray(np.meshgrid(spaced_values,spaced_values,indexing="ij"))
x_cord = intersection_mesh[0,:,:]
y_cord = intersection_mesh[1,:,:] 
#intersection_coordinates = np.rec.fromarrays([x_cord,y_cord],names='x,y')
intersection_coordinates = np.vstack(([x_cord.T], [y_cord.T])).T
peak_mask = img > 0
intensities = img[peak_mask]
intersection_coordinates = intersection_coordinates[peak_mask]

import pycrystem as pc

alpha = pc.diffraction_generator.DiffractionSimulation(coordinates=intersection_coordinates,intensities=intensities)
alpha = alpha.as_signal(144,0.08,2.5)
diffraction = pc.ElectronDiffraction(alpha.data)
diffraction.plot()
