import numpy as np
import scipy.sparse as sp
from scipy.signal import convolve2d


def smooothing(mat, h=1):
    sparse = isinstance(mat, sp.csr_matrix)
    if sparse:
        mat = mat.toarray()
    kernel = np.ones((h+1, h+1)) / (h+1)**2
    new_mat = convolve2d(mat, kernel, mode='same')
    if sparse:
        return sp.csr_matrix(new_mat)
    else:
        return new_mat


def random_walk():
    pass

