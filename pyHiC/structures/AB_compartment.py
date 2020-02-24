import numpy as np
from scipy.linalg import eigh
from ..normalization import normalization


def AB_compartment(mat, n_th_eigenvalue=1):
    assert n_th_eigenvalue in [1, 2]

    norm_mat = normalization(mat, 'OE')
    # Remove zeros
    sm = np.sum(norm_mat, axis=1)
    zeros = []
    for i in range(len(sm)):
        if sm[i] == 0:
            zeros.append(i)
    if len(zeros) != 0:
        norm_mat = np.delete(norm_mat, zeros, axis=0)
        norm_mat = np.delete(norm_mat, zeros, axis=1)

    degree_ = np.diag(1 / np.sqrt(np.sum(norm_mat, axis=1)))

    norm_laplacian = np.eye(len(norm_mat)) - degree_.dot(norm_mat).dot(degree_)
    eigen_val, eigen_vec = eigh(norm_laplacian, eigvals=(0, 1))

    ab_comp = eigen_vec[:, n_th_eigenvalue - 1]
    for zero in zeros:
        ab_comp = np.insert(ab_comp, zero, 0)

    return ab_comp

