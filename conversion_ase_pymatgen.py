import pymatgen as pmg
from pymatgen.io import ase

si = pmg.Element("Si")
lattice = pmg.Lattice.cubic(5.431)
silicon = pmg.Structure.from_spacegroup("Fd-3m",lattice, [si], [[0, 0, 0]])

silicon_ase = ase.AseAtomsAdaptor.get_atoms(silicon)
