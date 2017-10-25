import numpy as np
import hyperspy.api as hs
import pycrystem as pc
import pymatgen as pmg
from pymatgen.transformations.standard_transformations import DeformStructureTransformation
from pycrystem.indexation_generator import IndexationGenerator
from scipy.constants import pi

## Creating our structure
si = pmg.Element("Si")
lattice = pmg.Lattice.cubic(5.431)
structure = pmg.Structure.from_spacegroup("Fd-3m",lattice, [si], [[0, 0, 0]])
structure.make_supercell([3,3,3])

## Creating our electron diffraction set up
ediff = pc.ElectronDiffractionCalculator(300., 0.025)

diff_dat = ediff.calculate_ed_data(structure,algorithm='multi-slice',wave_size=144,num_slices=10)
dpi = diff_dat.as_signal(512, 0.08, 2.5)
diffraction = pc.ElectronDiffraction(dpi.data)
diffraction.plot()

