import pymatgen as pmg 
from pymatgen.io.cif import *

jif = CifParser('/home/phillip/Documents/code/phd/sandbox/azif.cif')

jif = jif.get_structures()