import numpy as np
from scipy.sparse import coo_matrix


def normalization(mat, method, **kwargs):
    if method not in ['log', 'power', 'OE', 'KR', 'VC', 'VC_SQRT']:
        print("Normalization operation not in ['OE', 'KR', 'VC', 'VC_SQRT']. Normalization omitted.")

    sparse = isinstance(mat, coo_matrix)
    if sparse:
        mat = mat.toarray()

    if method == 'log':
        if 'base' not in kwargs:
            base = np.e
        else:
            base = kwargs['base']
        mat = np.log(mat + 1) / np.log(base)

    if method == 'power':
        if 'power' not in kwargs:
            p = 0.5
        else:
            p = kwargs['power']
        mat = np.power(mat, p)

    if method == 'OE':
        averages = np.array([np.mean(mat[i:, :len(mat)-i]) for i in range(len(mat))])
        averages = np.where(averages == 0, 1, averages)
        for i in range(len(mat)):
            for j in range(len(mat)):
                d = abs(i - j)
                mat[i, j] = mat[i, j] / averages[d]

    if method == 'VC':
        sm = np.sum(mat, axis=0)
        sm = np.where(sm == 0, 1, sm)
        sm_v = np.tile(sm, (len(sm), 1))
        sm_c = sm_v.T
        mat = mat / sm_c / sm_v

    if method == 'VC_SQRT':
        sm = np.sum(mat, axis=0)
        sm = np.where(sm == 0, 1, sm)
        sm = np.sqrt(sm)
        sm_v = np.tile(sm, (len(sm), 1))
        sm_c = sm_v.T
        mat = mat / sm_c / sm_v

    if method == 'KR':
        bias = np.mean(mat) / 1e3
        # Remove all-zero rows and columns
        sm = np.sum(mat, axis=0)
        zeros = []
        for i in range(len(sm)):
            if sm[i] == 0:
                zeros.append(i)
        mat = np.delete(mat, zeros, axis=0)
        mat = np.delete(mat, zeros, axis=1)

        # Iteration
        x = np.random.random(size=len(mat))
        k = 0
        while True:
            k += 1
            aa = np.diag(x).dot(mat) + np.diag(mat.dot(x))
            aa = np.linalg.inv(aa)
            bb = np.diag(x).dot(mat).dot(x) - np.ones(x.shape)
            delta = aa.dot(bb)
            new_x = x - aa.dot(bb)

            max_error = np.max(np.abs(delta))
            # print(f'Iteration: {k}, Max Error: {max_error}')
            if max_error < bias:
                break
            else:
                x = new_x

        # Normalization
        dg = np.diag(new_x)
        mat = dg.dot(mat).dot(dg)

        # Put all-zero rows and columns back
        for zero in zeros:
            mat = np.insert(mat, zero, 0, axis=0)
            mat = np.insert(mat, zero, 0, axis=1)

    if sparse:
        mat = coo_matrix(mat)
    return mat

