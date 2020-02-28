import sys
sys.path.insert(0, '../')
print(sys.path)
import numpy as np
from pyHiC.loading import load_HiC
from pyHiC.structures import AB_compartment
from pyHiC.normalization import normalization
from pyHiC.visualization import visualize_HiC_triangle
from pyHiC.reproducibility import HiCRep



x = load_HiC('/nfs/turbo/umms-drjieliu/proj/4dn/data/bulkHiC/H1-hESC/processed/chr16_1kb.txt',
             format='short', chromosome='chr7', resolution=5000, sparse=False,
             start_pos=22000000, end_pos=23500000
             )
x = np.log(x + 1)
# np.save('chr16_25M_40M_5kb.npy', x)

# x = np.load('chr16_25M_40M_5kb.npy')
print(x.shape)
visualize_HiC_triangle(x, 'chr16_17M_32M_5kb.png', vmax=4, colorbar=False)

# y = x[900:1200, 900:1200]
# visualize_HiC_triangle(y, 'chr7_29.5M_31M_5kb.png', vmax=5, colorbar=False)

# z = x[995:1045, 995:1045]
# visualize_HiC_triangle(z, 'chr7_29.975M_30.225M_5kb.png', vmax=5, colorbar=False)

print(HiCRep(x, y))

