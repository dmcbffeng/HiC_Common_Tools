import numpy as np
import scipy.sparse as sp
from ..utils import smooothing


def HiCRep(HiC1, HiC2, n_strata=10, h=1):
    """
    Calculate the similarity of two contact maps with HiCRep score
    Args:
        HiC1, HiC2 (numpy.array or sp.csr_matrix): two Hi-C contact maps
        n_strata (int): Use first n strata (closest to the diagonal)
        h (int): size of smoothing window

    Return:
         HiCRep score (float)
    """
    assert HiC1.shape == HiC2.shape
    sparse1, sparse2 = isinstance(HiC1, sp.csr_matrix), isinstance(HiC2, sp.csr_matrix)
    if h != 0:
        HiC1, HiC2 = smooothing(HiC1, h), smooothing(HiC2, h)

    def _strata(_mat, _n_strata, _sparse):
        _l = _mat.shape[0]
        _str = [np.zeros((_l - i,)) for i in range(_n_strata)]
        if _sparse:
            _mat = _mat.tocoo()
            for r, c, v in zip(_mat.row, _mat.col, _mat.data):
                if abs(r - c) < _n_strata and r <= c:
                    _str[abs(r - c)][r] = v
        else:
            for i in range(_n_strata):
                _str[i] = np.diag(_mat[i:, :_l-i])
        return _str

    strata1, strata2 = _strata(HiC1, n_strata, sparse1), _strata(HiC2, n_strata, sparse2)
    print(strata1[0], strata2[0])
    std1, std2 = np.array([np.std(_s1) for _s1 in strata1]), np.array([np.std(_s2) for _s2 in strata2])
    lengths = np.arange(HiC1.shape[0], HiC1.shape[0] - n_strata, -1)
    weights = std1 * std2 * lengths
    Pearson_corrs = np.array([np.corrcoef(_s1, _s2)[0][1] for _s1, _s2 in zip(strata1, strata2)])
    score = np.sum(weights * Pearson_corrs / np.sum(weights))
    return score

