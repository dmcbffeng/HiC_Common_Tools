"""
Utils for loading contact maps in different data formats.
Used by other functions.
"""

import numpy as np
from scipy.sparse import coo_matrix, load_npz


def file_line_generator(file, chrom=None, header=False, format=None):
    if len(format) == 3:
        format = [0, format[0], 0, format[1], format[2]]
    count = 0
    f = open(file)
    if header:
        next(f)
    for line in f:
        count += 1
        if count % 100000 == 0:
            print('Line: ', count)

        lst = line.strip().split()
        if len(format) not in [4, 5]:
            raise ValueError('Wrong custom format!')
        if format[0] != 0 and format[2] != 0:
            c1, c2 = lst[format[0]-1], lst[format[2]-1]
            if c1 != chrom or c2 != chrom:
                continue
        p1, p2 = int(lst[format[1]-1]), int(lst[format[3]-1])
        if len(format) == 4:
            v = 1
        else:
            v = float(lst[format[4]-1])
        yield p1, p2, v


def load_HiC(file, format=None, custom_format=None, header=False,
             chromosome=None, start_pos=0, end_pos=-1,
             resolution=10000, sparse=False):
    if format in ['npy', 'npz']:
        if format == 'npy':
            mat = np.load(file)
        else:
            mat = load_npz(file)
            mat = mat.toarray()

    else:
        if format in ['short', 'Short']:
            gen = file_line_generator(file, format=[0, 1, 0, 2, 3])
        elif format in ['long', 'Long']:
            gen = file_line_generator(file, chrom=chromosome, format=[1, 2, 3, 4, 5])
        elif format in ['NoScore', 'noscore']:
            gen = file_line_generator(file, chrom=chromosome, format=[1, 2, 3, 4, 0])
        elif format is None:
            if custom_format is None:
                raise ValueError('Please provide file format!')
            else:
                if isinstance(custom_format, int):
                    custom_format = [int(elm) for elm in str(custom_format)]
                gen = file_line_generator(file, header=header, chrom=chromosome, format=custom_format)
        else:
            raise ValueError('Unrecognized format: ' + format)

        if end_pos != -1:
            size = int(np.ceil((end_pos - start_pos) / resolution))
        else:
            size = 100  # Current size set to be 100
            # But record max position all the time, double the size when necessary

        max_pos = 0
        mat = np.zeros((size, size))
        for p1, p2, val in gen:
            if end_pos != -1:
                if p1 >= end_pos or p2 >= end_pos:
                    continue
            p1, p2 = (p1 - start_pos) // resolution, (p2 - start_pos) // resolution
            if end_pos == -1:
                max_pos = max(p1, p2, max_pos)
                while max_pos >= size:
                    size *= 2
                    print('New Size: ', size)
                    new_mat = np.zeros((size, size))
                    new_mat[:len(mat), :len(mat)] = mat
                    mat = new_mat

            mat[p1, p2] += val
            if p1 != p2:
                mat[p2, p1] += val
        if end_pos == -1:
            mat = mat[:max_pos+1, :max_pos+1]

    if sparse:
        mat = coo_matrix(mat)
    return mat

