"""
Visualize a given region of HiC contact map or compare two contact maps.
"""
import numpy as np
import scipy.sparse as sp
import matplotlib.pyplot as plt
import seaborn as sns


def visualize_HiC_triangle(HiC, output, fig_size=(12, 6.5),
                           vmin=0, vmax=None, cmap='Reds', colorbar=True,
                           colorbar_orientation='vertical',
                           x_ticks=None, fontsize=24):
    """
        Visualize matched HiC and epigenetic signals in one figure
        Args:
            HiC (numpy.array): Hi-C contact map, only upper triangle is used.
            output (str): the output path. Must in a proper format (e.g., 'png', 'pdf', 'svg', ...).
            fig_size (tuple): (width, height). Default: (12, 8)
            vmin (float): min value of the colormap. Default: 0
            vmax (float): max value of the colormap. Will use the max value in Hi-C data if not specified.
            cmap (str or plt.cm): which colormap to use. Default: 'Reds'
            colorbar (bool): whether to add colorbar for the heatmap. Default: True
            colorbar_orientation (str): "horizontal" or "vertical". Default: "vertical"
            x_ticks (list): a list of strings. Will be added at the bottom. THE FIRST TICK WILL BE AT THE START OF THE SIGNAL, THE LAST TICK WILL BE AT THE END.
            fontsize (int): font size. Default: 24

        No return. Save a figure only.
        """
    if isinstance(HiC, sp.csr_matrix):
        HiC = HiC.toarray()

    N = len(HiC)
    coordinate = np.array([[[(x + y) / 2, y - x] for y in range(N + 1)] for x in range(N + 1)])
    X, Y = coordinate[:, :, 0], coordinate[:, :, 1]
    vmax = vmax if vmax is not None else np.max(HiC)

    fig, ax = plt.subplots(figsize=fig_size)
    im = plt.pcolormesh(X, Y, HiC, vmin=vmin, vmax=vmax, cmap=cmap)
    # plt.axis('off')
    plt.yticks([], [])
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    if x_ticks:
        tick_pos = np.linspace(0, N, len(x_ticks))
        ax.set_xticks(tick_pos)
        ax.set_xticklabels(x_ticks, fontsize=fontsize)
    else:
        ax.spines['bottom'].set_visible(False)

    plt.ylim([0, N])
    plt.xlim([0, N])

    if colorbar:
        if colorbar_orientation == 'horizontal':
            _left, _width, _bottom, _height = 0.7, 0.25, 0.75, 0.03
        elif colorbar_orientation == 'vertical':
            _left, _width, _bottom, _height = 0.9, 0.02, 0.3, 0.5
        else:
            raise ValueError('Wrong orientation!')
        cbar = plt.colorbar(im, cax=fig.add_axes([_left, _bottom, _width, _height]),
                            orientation=colorbar_orientation)
        cbar.ax.tick_params(labelsize=fontsize)
        cbar.outline.set_visible(False)

    plt.savefig(output)
    # plt.show()


def visualize_HiC_square(HiC, output, fig_size=(12, 6.5),
                         vmin=0, vmax=None, cmap='Reds', colorbar=True,
                         x_ticks=None, fontsize=24):
    """
        Visualize matched HiC and epigenetic signals in one figure
        Args:
            HiC (numpy.array): Hi-C contact map, only upper triangle is used.
            output (str): the output path. Must in a proper format (e.g., 'png', 'pdf', 'svg', ...).
            fig_size (tuple): (width, height). Default: (12, 8)
            vmin (float): min value of the colormap. Default: 0
            vmax (float): max value of the colormap. Will use the max value in Hi-C data if not specified.
            cmap (str or plt.cm): which colormap to use. Default: 'Reds'
            colorbar (bool): whether to add colorbar for the heatmap. Default: True
            x_ticks (list): a list of strings. Will be added at the bottom. THE FIRST TICK WILL BE AT THE START OF THE SIGNAL, THE LAST TICK WILL BE AT THE END.
            fontsize (int): font size. Default: 24

        No return. Save a figure only.
        """
    if isinstance(HiC, sp.csr_matrix):
        HiC = HiC.toarray()

    N = len(HiC)
    vmax = vmax if vmax is not None else np.max(HiC)

    plt.subplots(figsize=fig_size)
    if x_ticks:
        g = sns.heatmap(HiC, vmin=0, vmax=vmax, cmap='Reds', square=True,
                        xticklabels=x_ticks, yticklabels=x_ticks, cbar=colorbar)
        tick_pos = np.linspace(0, N, len(x_ticks))
        g.set_xticks(tick_pos)
        g.set_xticklabels(g.get_xmajorticklabels(), fontsize=fontsize, rotation=0)
        g.set_yticks(tick_pos)
        g.set_yticklabels(g.get_ymajorticklabels(), fontsize=fontsize)
    else:
        g = sns.heatmap(HiC, vmin=0, vmax=vmax, cmap='Reds', square=True, cbar=colorbar)
        g.set_xticks([])
        g.set_yticks([])
    plt.savefig(output)
    # plt.show()


if __name__ == '__main__':
    np.random.seed(1)
    hic = np.random.random((100, 100)) + 1
    hic = sp.csr_matrix(hic)
    visualize_HiC_square(hic, 'test.png', x_ticks=['1', '3', '5'], colorbar=False)


