"""
Utils for loading contact maps in different data formats.
Used by other functions.
"""

import numpy as np
import scipy.sparse as sp
import networkx as nx
import gzip as gz


def file_line_generator(file, chrom=None, header=False, format=None, gzip=False):
    if len(format) == 3:
        format = [0, format[0], 0, format[1], format[2]]
    count = 0
    f = gz.open(file) if gzip else open(file)
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
             resolution=10000, gzip=False, sparse=True):
    """
    Load the contact matrix of one chromosome (or part of one chromosome) from a HiC file

    Args:
        file: (str) file name;
        format: (str or None) default: None. "short" / "Short", "long" / "Long", "noscore" / "NoScore", "npy" or "npz". If customized, leave it "None".
        custom_format: (str or list or None) default: None. For customized input, provide the indices like "2356"".
        header: (bool or None) default: None. For customized input, whether the file has a header line.
        chromosome: (str) default: None. For formats other than "short", give the chromosome you would like to extract, eg. "chr1".
        start_pos & end_pos: (int) default: 0 and -1. (0: start, -1: end).
        resolution: (int) default: 10000.
        gzip (bool): whether zipped file. Default: False
        sparse: (bool) default: True. If True, store with scipy.sparse.csr_matrix; if false, with numpy.array.

    Return:
         HiC contact matrix (numpy.array or scipy.sparse.csr_matrix)
    """

    if format in ['npy', 'npz']:
        if format == 'npy':
            mat = np.load(file)
            if sparse:
                mat = sp.csr_matrix(mat)
        else:
            mat = sp.load_npz(file)
            if not sparse:
                mat = mat.toarray()
            if sparse and not isinstance(mat, sp.csr_matrix):
                mat = mat.tocsr()

    else:
        if format in ['short', 'Short']:
            gen = file_line_generator(file, format=[0, 1, 0, 2, 3], gzip=gzip)
        elif format in ['long', 'Long']:
            gen = file_line_generator(file, chrom=chromosome, format=[1, 2, 3, 4, 5], gzip=gzip)
        elif format in ['NoScore', 'noscore']:
            gen = file_line_generator(file, chrom=chromosome, format=[1, 2, 3, 4, 0], gzip=gzip)
        elif format is None:
            if custom_format is None:
                raise ValueError('Please provide file format!')
            else:
                if isinstance(custom_format, int):
                    custom_format = [int(elm) for elm in str(custom_format)]
                gen = file_line_generator(file, header=header, chrom=chromosome, format=custom_format, gzip=gzip)
        else:
            raise ValueError('Unrecognized format: ' + format)

        assert end_pos == -1 or end_pos > start_pos
        size = int(np.ceil((end_pos - start_pos) / resolution)) if end_pos != -1 else 1

        G = nx.Graph()
        for i in range(size):
            G.add_node(i)
            G.add_edge(i, i, weight=np.finfo(float).eps)

        for p1, p2, val in gen:
            if p1 < start_pos or p2 < start_pos:
                continue
            if end_pos != -1:
                if p1 >= end_pos or p2 >= end_pos:
                    continue
            p1, p2 = (p1 - start_pos) // resolution, (p2 - start_pos) // resolution
            if end_pos == -1:
                max_pos = max(p1, p2)
                if max_pos > size:
                    for i in range(size, max_pos + 1):
                        G.add_node(i)
                        G.add_edge(i, i, weight=np.finfo(float).eps)
                    size = max_pos + 1

            if G.has_edge(p1, p2):
                G[p1][p2]['weight'] += val
            else:
                G.add_edge(p1, p2, weight=val)

        mat = nx.adj_matrix(G)
        if not sparse:
            mat = mat.toarray()
    return mat

