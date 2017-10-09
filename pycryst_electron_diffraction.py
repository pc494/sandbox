
# coding: utf-8



#get_ipython().magic('matplotlib tk')

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

## Creating our electron diffraction set up
ediff = pc.ElectronDiffractionCalculator(300., 0.025)

#help(ediff.calculate_ed_data)

diff_dat = ediff.calculate_ed_data(structure,algorithm='multi-slice',reciprocal_radius=2.5)

dpi = diff_dat.as_signal(512, 0.02, 2.5)

diffraction = pc.ElectronDiffraction(dpi.data)

diffraction.plot()

