"""
Visualize matched HiC and epigenetic signals in one figure

@author: Fan Feng
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec


def visualize_HiC_epigenetics(HiC, epis, output, fig_width=12.0,
                              vmin=0, vmax=None, cmap='Reds', colorbar=True,
                              colorbar_orientation='vertical',
                              epi_labels=None, x_ticks=None, fontsize=24,
                              epi_colors=None, epi_yaxis=True,
                              heatmap_ratio=0.6, epi_ratio=0.1,
                              interval_after_heatmap=0.05, interval_between_epi=0.01,):
    """
    Visualize matched HiC and epigenetic signals in one figure
    Args:
        HiC (numpy.array): Hi-C contact map, only upper triangle is used.
        epis (list): epigenetic signals
        output (str): the output path. Must in a proper format (e.g., 'png', 'pdf', 'svg', ...).
        fig_width (float): the width of the figure. Then the height will be automatically calculated. Default: 12.0
        vmin (float): min value of the colormap. Default: 0
        vmax (float): max value of the colormap. Will use the max value in Hi-C data if not specified.
        cmap (str or plt.cm): which colormap to use. Default: 'Reds'
        colorbar (bool): whether to add colorbar for the heatmap. Default: True
        colorbar_orientation (str): "horizontal" or "vertical". Default: "vertical"
        epi_labels (list): the names of epigenetic marks. If None, there will be no labels at y axis.
        x_ticks (list): a list of strings. Will be added at the bottom. THE FIRST TICK WILL BE AT THE START OF THE SIGNAL, THE LAST TICK WILL BE AT THE END.
        fontsize (int): font size. Default: 24
        epi_colors (list): colors of epigenetic signals
        epi_yaxis (bool): whether add y-axis to epigenetic signals. Default: True
        heatmap_ratio (float): the ratio of (heatmap height) and (figure width). Default: 0.6
        epi_ratio (float): the ratio of (1D epi signal height) and (figure width). Default: 0.1
        interval_after_heatmap (float): the ratio of (interval between heatmap and 1D signals) and (figure width). Default: 0.05
        interval_between_epi (float): the ratio of (interval between 1D signals) and (figure width). Default: 0.01

    No return. Save a figure only.
    """
    # Make sure the lengths match
    len_epis = [len(epi) for epi in epis]
    if max(len_epis) != min(len_epis) or max(len_epis) != len(HiC):
        raise ValueError('Size not matched!')
    N = len(HiC)

    # Define the space for each row (heatmap - interval - signal - interval - signal ...)
    rs = [heatmap_ratio, interval_after_heatmap] + [epi_ratio, interval_between_epi] * len(epis)
    rs = np.array(rs[:-1])
    rs = rs / np.sum(rs)  # normalize to 1 (ratios)

    # Calculate figure height
    fig_height = fig_width * np.sum(rs)
    fig = plt.figure(figsize=(fig_width, fig_height))

    # Split the figure into rows with different heights
    gs = GridSpec(len(rs), 1, height_ratios=rs)

    # Ready for plotting heatmap
    ax0 = plt.subplot(gs[0, :])
    # Define the rotated axes and coordinates
    coordinate = np.array([[[(x + y) / 2, y - x] for y in range(N + 1)] for x in range(N + 1)])
    X, Y = coordinate[:, :, 0], coordinate[:, :, 1]
    # Plot the heatmap
    vmax = vmax if vmax is not None else np.max(HiC)
    im = ax0.pcolormesh(X, Y, HiC, vmin=vmin, vmax=vmax, cmap=cmap)
    ax0.axis('off')
    ax0.set_ylim([0, N])
    ax0.set_xlim([0, N])
    if colorbar:
        if colorbar_orientation == 'horizontal':
            _left, _width, _bottom, _height = 0.12, 0.25, 1 - rs[0] * 0.25, rs[0] * 0.03
        elif colorbar_orientation == 'vertical':
            _left, _width, _bottom, _height = 0.9, 0.02, 1 - rs[0] * 0.7, rs[0] * 0.5
        else:
            raise ValueError('Wrong orientation!')
        cbar = plt.colorbar(im, cax=fig.add_axes([_left, _bottom, _width, _height]),
                            orientation=colorbar_orientation)
        cbar.ax.tick_params(labelsize=fontsize)
        cbar.outline.set_visible(False)

    # print(rs/np.sum(rs))
    # Ready for plotting 1D signals
    if epi_labels:
        assert len(epis) == len(epi_labels)
    if epi_colors:
        assert len(epis) == len(epi_colors)

    for i, epi in enumerate(epis):
        # print(epi.shape)
        ax1 = plt.subplot(gs[2 + 2 * i, :])
        if epi_colors:
            ax1.fill_between(np.arange(N), 0, epi, color=epi_colors[i])
        else:
            ax1.fill_between(np.arange(N), 0, epi)
        if not epi_yaxis:
            ax1.set_yticks([])
            ax1.set_yticklabels([])
            ax1.spines['left'].set_visible(False)
        else:
            ax1.tick_params(labelsize=fontsize)
        if i != len(epis) - 1:
            ax1.set_xticks([])
            ax1.set_xticklabels([])
        # ax1.axis('off')
        # ax1.xaxis.set_visible(True)
        # plt.setp(ax1.spines.values(), visible=False)
        # ax1.yaxis.set_visible(True)
        ax1.spines['right'].set_visible(False)
        ax1.spines['top'].set_visible(False)
        ax1.spines['bottom'].set_visible(False)
        ax1.set_xlim([-0.5, N-0.5])
        if epi_labels:
            ax1.set_ylabel(epi_labels[i], fontsize=fontsize)
    ax1.spines['bottom'].set_visible(True)
    if x_ticks:
        tick_pos = np.linspace(0, N - 1, len(x_ticks))   # 这个坐标其实是不对的 差1个bin 但是为了ticks好看只有先这样了
        ax1.set_xticks(tick_pos)
        ax1.set_xticklabels(x_ticks, fontsize=fontsize)
    else:
        ax1.set_xticks([])
        ax1.set_xticklabels([])

    plt.savefig(output)


if __name__ == '__main__':
    np.random.seed(1)
    hic = np.random.random((100, 100)) + 1
    s1, s2 = np.random.random((100,)) + 1, np.random.random((100,)) + 1
    visualize_HiC_epigenetics(hic, [s1, s2], 'test.png', interval_after_heatmap=0.,
                              interval_between_epi=0., epi_ratio=0.15,
                              vmax=2, colorbar_orientation='vertical',
                              epi_colors=['red', 'green'], x_ticks=['1.1Mb', '1.2Mb', '1.3Mb'], epi_yaxis=False
                              )


