import numpy as np
import scipy
import scipy.sparse as sp
from scipy.linalg import inv


x = sp.csr_matrix(([1,2,3], ([0,1,2], [0,1,2])), shape=(3,3))
print(x)