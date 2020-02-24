import sys
sys.path.insert(0, '../')
print(sys.path)
import numpy as np
from pyHiC.loading import load_HiC
from pyHiC.structures import AB_compartment
from pyHiC.normalization import normalization
from pyHiC.visualization import visualize_HiC_triangle


x = load_HiC('hESC_chr1_100kb.txt', format='short', chromosome='chr1', resolution=1000000, sparse=False)
x = np.log(x + 1)


visualize_HiC_triangle(x, 'test.png', vmax=15)


