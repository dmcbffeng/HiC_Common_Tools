import numpy as np
import scipy.sparse as sp
import time
import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)


def normalization(mat, method, **kwargs):
    if method.lower() not in ['oe', 'kr', 'vc', 'vc_sqrt', 'ic']:
        print("Normalization operation not in ['OE', 'KR', 'VC', 'VC_SQRT', 'IC']. Normalization omitted.")

    sparse = isinstance(mat, sp.csr_matrix)

    if method.lower() == 'oe':
        if sparse:
            l = mat.shape[0]
            sm_strata = np.zeros((l,))
            new_mat = mat.tocoo()
            for r, c, v in zip(new_mat.row, new_mat.col, new_mat.data):
                sm_strata[abs(r - c)] += v
            mean_strata = sm_strata / np.arange(l, 0, -1)
            vs = []
            for r, c, v in zip(new_mat.row, new_mat.col, new_mat.data):
                vs.append(v / mean_strata[abs(r - c)])
            mat = sp.csr_matrix((vs, (new_mat.row, new_mat.col)), shape=mat.shape)

        else:
            averages = np.array([np.mean(np.diag(mat[i:, :len(mat)-i])) for i in range(len(mat))])
            averages = np.where(averages == 0, 1, averages)
            for i in range(len(mat)):
                for j in range(i, len(mat)):
                    val = mat[i, j] / averages[abs(i - j)]
                    mat[i, j] = val
                    mat[j, i] = val

    if method.lower() == 'vc':
        if sparse:
            sm = 1 / mat.sum(axis=0).toarray().flatten()
            sm_ = sp.diags(np.where(sm == 0, 1, sm), 0)
            mat = sm_ * mat * sm_
        else:
            sm = np.sum(mat, axis=0)
            sm_ = np.diag(1 / np.where(sm == 0, 1, sm))
            mat = sm_.dot(mat).dot(sm_)

    if method.lower() == 'vc_sqrt':
        if sparse:
            sm = 1 / np.sqrt(mat.sum(axis=0).toarray().flatten())
            sm_ = sp.diags(np.where(sm == 0, 1, sm), 0)
            mat = sm_ * mat * sm_
        else:
            sm = np.sqrt(np.sum(mat, axis=0))
            sm_ = np.diag(1 / np.where(sm == 0, 1, sm))
            mat = sm_.dot(mat).dot(sm_)

    if method.lower() == 'kr':
        if sparse:
            mat = mat.toarray()

        tolerance = kwargs.pop('tolerance', 1e-5)
        verbose = kwargs.pop('verbose', 1)
        max_iteration = kwargs.pop('max_iteration', 50)
        # Remove all-zero rows and columns
        sm = mat.sum(axis=0).toarray().flatten() if sparse else np.sum(mat, axis=0)
        zeros = []
        for i in range(len(sm)):
            if sm[i] == 0:
                zeros.append(i)

        mat = np.delete(mat, zeros, axis=0)
        mat = np.delete(mat, zeros, axis=1)

        # Iteration
        x = np.random.random(size=len(mat))
        new_x = None
        if verbose:
            print("[KR Norm] starting iterative correction")
        start_time = time.time()
        for _iter in range(max_iteration):
            aa = np.diag(x).dot(mat) + np.diag(mat.dot(x))
            aa = np.linalg.inv(aa)
            bb = np.diag(x).dot(mat).dot(x) - np.ones(x.shape)
            delta = aa.dot(bb)
            new_x = x - aa.dot(bb)

            max_error = np.sum(np.abs(delta))
            # print(f'Iteration: {k}, Max Error: {max_error}')
            if verbose:
                if _iter % 2 == 1:
                    end_time = time.time()
                    estimated = (float(max_iteration - _iter - 1) * (end_time - start_time)) / (_iter + 1)
                    m, sec = divmod(estimated, 60)
                    h, m = divmod(m, 60)
                    print("[KR Norm] pass {} Estimated time {:.0f}:{:.0f}:{:.0f}".format(_iter, h, m, sec))
                    # print("max delta - 1 = {} ".format(deviation))

            if max_error < tolerance:
                if verbose:
                    print("[KR Norm] {} iterations used\n".format(_iter + 1))
                break
            else:
                x = new_x
        else:
            print("[KR Norm] Max {} iterations reached. Iteration stopped.\n".format(max_iteration + 1))

        # Normalization
        dg = np.diag(new_x)
        mat = dg.dot(mat).dot(dg)

        # Put all-zero rows and columns back
        for zero in zeros:
            mat = np.insert(mat, zero, 0, axis=0)
            mat = np.insert(mat, zero, 0, axis=1)

        if sparse:
            mat = sp.csr_matrix(mat)

    if method.lower() == 'ic':
        mat = iterativeCorrection(mat, kwargs.pop('max_iteration', 50),
                                  kwargs.pop('tolerance=1e-5', kwargs.pop('verbose', 1)))

    return mat


def iterativeCorrection(matrix, max_iteration=50, tolerance=1e-5, verbose=1):
    """
    adapted from cytonised version in mirnylab
    original code from: ultracorrectSymmetricWithVector
    https://bitbucket.org/mirnylab/mirnylib/src/924bfdf5ed344df32743f4c03157b0ce49c675e6/mirnylib/numutils_new.pyx?at=default
    Main method for correcting DS and SS read data.
    Possibly excludes diagonal.
    By default does iterative correction, but can perform an M-time correction
    :param matrix: a scipy sparse matrix
    :param tolerance: Tolerance is the maximum allowed relative
                      deviation of the marginals.
    """
    sparse = isinstance(matrix, sp.csr_matrix)
    total_bias = np.ones(matrix.shape[0], 'float64')

    if np.isnan(np.sum(matrix)):
        print("[iterative correction] the matrix contains nans, they will be replaced by zeros.")
        matrix.data[np.isnan(matrix.data)] = 0

    matrix = matrix.astype(float)
    W = matrix.tocoo()

    if np.abs(matrix - matrix.T).mean() / (1. * np.abs(matrix.mean())) > 1e-10:
        raise ValueError("Please provide symmetric matrix!")

    start_time = time.time()
    if verbose:
        print("starting iterative correction")
    for iternum in range(max_iteration):
        iternum += 1
        s = np.array(W.sum(axis=1)).flatten()
        mask = (s == 0)
        s = s / np.mean(s[~mask])

        total_bias *= s
        deviation = np.abs(s - 1).max()

        s = 1.0 / s

        # The following code  is an optimization of this
        # for i in range(N):
        #     for j in range(N):
        #         W[i,j] = W[i,j] / (s[i] * s[j])

        W.data *= np.take(s, W.row)
        W.data *= np.take(s, W.col)
        if np.any(W.data > 1e100):
            print("*Error* matrix correction is producing extremely large values. "
                      "This is often caused by bins of low counts. Use a more stringent "
                      "filtering of bins.")
            exit(1)
        if verbose:
            if iternum % 5 == 0:
                end_time = time.time()
                estimated = (float(max_iteration - iternum - 1) * (end_time - start_time)) / (iternum + 1)
                m, sec = divmod(estimated, 60)
                h, m = divmod(m, 60)
                print("pass {} Estimated time {:.0f}:{:.0f}:{:.0f}".format(iternum + 1, h, m, sec))
                print("max delta - 1 = {} ".format(deviation))

        if deviation < tolerance:
            print("[iterative correction] {} iterations used\n".format(iternum + 1))
            break

    # scale the total bias such that the sum is 1.0
    corr = total_bias[total_bias != 0].mean()
    total_bias /= corr
    W.data = W.data * corr * corr
    if np.any(W.data > 1e10):
        print("*Error* matrix correction produced extremely large values. "
                  "This is often caused by bins of low counts. Use a more stringent "
                  "filtering of bins.")
        exit(1)

    if sparse:
        return W.tocsr()
    else:
        return W.toarray()


