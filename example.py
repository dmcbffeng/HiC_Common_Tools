from utils import load_HiC
# from visualize import *
from normalization import *
from AB_compartment import *


HiC_mat = load_HiC('data/rep1_chr1.txt', format='short', resolution=250000, sparse=False)
print(HiC_mat.shape)


