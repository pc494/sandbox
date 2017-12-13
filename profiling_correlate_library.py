import pycrystem as pc
import numpy as np
from hyperspy.signals import BaseSignal
from heapq import nlargest
from operator import itemgetter
from pycrystem.utils import correlate

dp = pc.load('/home/pc494/Documents/data/Jandrek/QD-03.blo')
dp = dp.inav[0:2,0:2]
dp.apply_affine_transformation(D = np.array([[0.99, 0.00, 0.00],
                                             [0.00, 0.69, 0.00],
                                             [0.00, 0.00, 1.00]]))
dp.set_calibration(0.032)
structure = pc.Structure.from_file("/home/pc494/Documents/data/Jandrek/CsPbBr3_1.cif")

rot_array = np.loadtxt('/home/pc494/Documents/data/Jandrek/mmm_grid_euler.bin')
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
### copy/paste job
@profile
def correlate_library(image, library, n_largest):
    """Correlates all simulated diffraction templates in a DiffractionLibrary
    with a particular experimental diffraction pattern (image) stored as a
    numpy array.
    """
    i=0
    out_arr = np.zeros((n_largest * len(library),5))
    for key in library.keys():
        if n_largest:
            pass
        else:
            n_largest=len(library[key])
        correlations = dict()
        for orientation, diffraction_pattern in library[key].items():
            correlation = correlate(image, diffraction_pattern)
            correlations[orientation] = correlation
        res = nlargest(n_largest, correlations.items(), key=itemgetter(1))
        for j in np.arange(n_largest):
            out_arr[j + i*n_largest][0] = i
            out_arr[j + i*n_largest][1] = res[j][0][0]
            out_arr[j + i*n_largest][2] = res[j][0][1]
            out_arr[j + i*n_largest][3] = res[j][0][2]
            out_arr[j + i*n_largest][4] = res[j][1]
        i = i + 1
    return out_arr


class MatchingResults(BaseSignal):
    _signal_type = "matching_results"
    _signal_dimension = 2

    def __init__(self, *args, **kwargs):
        BaseSignal.__init__(self, *args, **kwargs)
        self.axes_manager.set_signal_dimension(2)

class IndexationGenerator():
    """Generates an indexer for data using a number of methods.
    Parameters
    ----------
    signal : ElectronDiffraction
        The signal of electron diffraction patterns to be indexed.
    library : DiffractionLibrary
        The library of simulated diffraction patterns for indexation
    """
    def __init__(self, signal, library):
        self.signal = signal
        self.library = library

    def correlate(self,
                  n_largest=5,
                  *args, **kwargs):
        """Correlates the library of simulated diffraction patterns with the
        electron diffraction signal.
        Parameters
        ----------
        n_largest : integer
            The n orientations with the highest correlation values are returned.
        *args/**kwargs : keyword arguments
            Keyword arguments passed to the HyperSpy map() function. Important
            options include...
        Returns
        -------
        matching_results : ndarray
            Numpy array with the same shape as the the navigation axes of the
            electron diffraction signal containing correlation results for each
            diffraction pattern.
        """
        signal = self.signal
        library = self.library
        matching_results = signal.map(correlate_library,
                                      library=library,
                                      n_largest=n_largest,
                                      inplace=False,
                                      *args, **kwargs)
        return MatchingResults(matching_results)

indexer = IndexationGenerator(dp, library)
match_results = indexer.correlate(parallel = True,n_largest=3)
#print(match_results.inav[0,0].data)