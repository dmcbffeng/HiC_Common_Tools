"""
Visualize a given region of HiC contact map or compare two contact maps.
"""
import numpy as np
from scipy.sparse import coo_matrix
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import seaborn as sns


def visualize_one_contact_map(mat, vmax, save_path=None):
    sparse = isinstance(mat, coo_matrix)
    if sparse:
        mat = mat.toarray()

    plt.figure(figsize=(9, 9))
    sns.heatmap(mat, vmin=0, vmax=vmax, square=True, cmap='Reds')
    if save_path:
        plt.savefig(save_path)
    plt.show()


def visualize_one_contact_map_with_coordinate(mat, vmax, st_pos, ed_pos, save_path=None):
    sparse = isinstance(mat, coo_matrix)
    if sparse:
        mat = mat.toarray()

    plt.figure(figsize=(10, 10))
    xtick = np.linspace(0, len(mat), 11)
    xticklabel = [str(pos / 1000000)[:5] for pos in np.linspace(st_pos, ed_pos, 11)]

    g = sns.heatmap(mat, vmin=0, vmax=vmax, cmap='Reds', square=True, yticklabels=xticklabel,
                    xticklabels=xticklabel)
    g.set_xticks(xtick)
    g.set_yticks(xtick)
    g.set_yticklabels(g.get_ymajorticklabels(), fontsize=16)
    g.set_xticklabels(g.get_xmajorticklabels(), fontsize=16)
    plt.xlabel('Genomic Coordinate (Mb)', fontsize=16)
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path)
    plt.show()


def visualize_two_contact_maps(mat1, mat2, vmax, save_path=None):
    sparse = isinstance(mat1, coo_matrix)
    if sparse:
        mat1 = mat1.toarray()

    sparse = isinstance(mat2, coo_matrix)
    if sparse:
        mat2 = mat2.toarray()

    ui = np.triu_indices(len(mat1), len(mat1)-1)
    mat2[ui] = mat1[ui]
    plt.figure(figsize=(9, 9))
    sns.heatmap(mat2, vmin=0, vmax=vmax, square=True, cmap='Reds')
    if save_path:
        plt.savefig(save_path)
    plt.show()




