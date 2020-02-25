import sys
sys.path.insert(0, '../')
print(sys.path)
import numpy as np
from pyHiC.loading import load_HiC
from pyHiC.structures import AB_compartment
from pyHiC.normalization import normalization
from pyHiC.visualization import visualize_HiC_triangle


x = load_HiC('D:\google_cloud\pancreas\chr7_1kb.txt',
             format='short', chromosome='chr7', resolution=10000,
             sparse=False, start_pos=0, end_pos=10000000)
x = np.log(x + 1)
np.save('chr7_first_10mb_10kb.npy', x)

visualize_HiC_triangle(x, 'test.png', vmax=15)


