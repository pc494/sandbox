import pyprismatic as pr
from pycrystem.utils.pyprismatic_io_utils import import_pyprismatic_data
from matplotlib import pyplot as plt
import numpy as np

threeD = import_pyprismatic_data()
twoD = np.squeeze(np.sum(threeD,axis=2))
plt.imshow(twoD)
plt.show() 
