# Commonly-used HiC Operations for Liu Lab (Python 3) #
## -- [NOT FINISHED YET] -- ##
@Author: Fan Feng

[1 Installation](#installation)\
[2 Loading Hi-C Contact Maps](#load-hi-c-contact-maps)\
[3 Normalization](#normalization)\
[4 Visualization](#visualization)\
[5 Structure Calling](#structure-calling)\
[6 Comparing Contact Maps](#comparing-contact-maps)


# Installation
**Required Packages**
- numpy
- scipy
- networkx
- matplotlib
- seaborn

**Install from GitHub**\
You can install the package with following command:
  ```console
    $ git clone https://github.com/dmcbffeng/HiC_Common_Tools.git
    $ cd scHiCTools
    $ python setup.py install
  ```


# Load Hi-C Contact Maps
 **Supported Formats**
 - npy: numpy.array / numpy.matrix
 - npz: scipy.sparse.coo_matrix / csr_matrix
 - Short
 ```
 <position1> <position2> <score>
 ```
 - Long
 ```
 <chromosome1> <position1> <chromosome2> <position2> <score>
 ```
 - NoScore
 ```
 <chromosome1> <position1> <chromosome2> <position2>
 ```
 - .hic format: JuiceTools can process it into "short" format with "dump" function.
 - Other formats: Simply give the indices (start from 1) in the order of
 "chromosome1 - position1 - chromosome2 - position2 - score" or
 "chromosome1 - position1 - chromosome2 - position2" or
 "position1 - position2 - score".
 For example, you can provide "2356" or [2, 3, 5, 6] if the file takes this format:
 ```
 <name> <chromosome1> <position1> <frag1> <chromosome2> <position2> <frag2> <strand1> <strand2>
 contact_1 chr1 3000000 1 chr1 3001000 1 + -
 ```
 
 **Load into (Sparse) Matrices**
 ```console
 >>> from pyHiC.loading import load_HiC
 >>> HiC_mat = load_HiC(
 ...         file='ESC_chr1.txt', format='short',
 ...         custom_format=None, header=False,
 ...         chromosome=None, start_pos=0, end_pos=-1,
 ...         resolution=500000, sparse=True)
 ```
 - file: (str) file name;
 - format: (str or None) default: None. "short" / "Short", "long" / "Long", "noscore" / "NoScore", "npy" or "npz". If customized, leave it "None". 
 - custom_format: (str or list or None) default: None. For customized input, provide the indices like "2356"".
 - header: (bool or None) default: None. For customized input, whether the file has a header line.
 - chromosome: (str) default: None. For formats other than "short", give the chromosome you would like to extract, eg. "chr1".
 - start_pos & end_pos: (int) default: 0 and -1. (0: start, -1: end).
 - resolution: (int) default: 10000.
 - gzip (bool): whether zipped file. Default: False
 - sparse: (bool) default: True. If True, store with scipy.sparse.csr_matrix; if false, with numpy.array.
 

# Normalization
 ```config
 >>> from pyHiC.normalization import normalization
 >>> normalized_mat_1 = normalization(HiC_mat, method='log', base=10)
 >>> normalized_mat_2 = normalization(HiC_mat, method='VC_SQRT')
 ```
 Normalize Hi-C contact maps, return normalized map.
 - method (str):
   - "OE": each value divided by the average of its corresponding strata (diagonal line)
   - "VC": each value divided by the sum of corresponding row then
   divided by the sum of corresponding column
   - "VC_SQRT": each value divided by the sqrt of the sum of corresponding row then
   divided by the sqrt of the sum of corresponding column
   - "KR": the sum of each row / column is one
   - "IC": iterative correction
   
 For KR and IC normalization, optional arguments include:
   - max_iteration (int): default: 50
   - tolerance (float): default: 1e-5
   - verbose (int, 1 or 0): whether print iteration information. default: 1 

# Visualization
 ```config
 >>> from pyHiC.visualization import *
 >>> visualize_HiC_epigenetics(HiC_mat, epis, 'output0.png', fig_width=12.0,
 ...        vmin=0, vmax=None, cmap='Reds', colorbar=True,
 ...        colorbar_orientation='vertical',
 ...        epi_labels=None, x_ticks=None, fontsize=24,
 ...        epi_colors=None, epi_yaxis=True,
 ...        heatmap_ratio=0.6, epi_ratio=0.1,
 ...        interval_after_heatmap=0.05, interval_between_epi=0.01,)
 ```
 Visualize matched HiC and epigenetic signals in one figure.
 Then save the figure as a file.
 - HiC (numpy.array): Hi-C contact map, only upper triangle is used.
 - epis (list): epigenetic signals
 - output (str): the output path. Must in a proper format (e.g., 'png', 'pdf', 'svg', ...).
 - fig_width (float): the width of the figure. Then the height will be automatically calculated. Default: 12.0
 - vmin (float): min value of the colormap. Default: 0
 - vmax (float): max value of the colormap. Will use the max value in Hi-C data if not specified.
 - cmap (str or plt.cm): which colormap to use. Default: 'Reds'
 - colorbar (bool): whether to add colorbar for the heatmap. Default: True
 - colorbar_orientation (str): "horizontal" or "vertical". Default: "vertical"
 - epi_labels (list): the names of epigenetic marks. If None, there will be no labels at y axis.
 - x_ticks (list): a list of strings. Will be added at the bottom. THE FIRST TICK WILL BE AT THE START OF THE SIGNAL, THE LAST TICK WILL BE AT THE END.
 - fontsize (int): font size. Default: 24
 - epi_colors (list): colors of epigenetic signals
 - epi_yaxis (bool): whether add y-axis to epigenetic signals. Default: True
 - heatmap_ratio (float): the ratio of (heatmap height) and (figure width). Default: 0.6
 - epi_ratio (float): the ratio of (1D epi signal height) and (figure width). Default: 0.1
 - interval_after_heatmap (float): the ratio of (interval between heatmap and 1D signals) and (figure width). Default: 0.05
 - interval_between_epi (float): the ratio of (interval between 1D signals) and (figure width). Default: 0.01

 ```config
 >>> from pyHiC.visualization import *
 >>> visualize_HiC_triangle(HiC, 'output1.png', fig_size=(12, 6.5),
 ...    vmin=0, vmax=None, cmap='Reds', colorbar=True,
 ...    colorbar_orientation='vertical',
 ...    x_ticks=None, fontsize=24)
 >>> visualize_HiC_triangle(HiC, 'output2.png', fig_size=(12, 6.5),
 ...    vmin=0, vmax=None, cmap='Reds', colorbar=True,
 ...    colorbar_orientation='vertical',
 ...    x_ticks=None, fontsize=24)
 ```
 Visualize one HiC contact map in triangle or square shape
  - HiC (numpy.array): Hi-C contact map, only upper triangle is used.
  - output (str): the output path. Must in a proper format (e.g., 'png', 'pdf', 'svg', ...).
  - fig_size (tuple): (width, height). Default: (12, 6.5) for triangle and (12, 12) for square
  - vmin (float): min value of the colormap. Default: 0
  - vmax (float): max value of the colormap. Will use the max value in Hi-C data if not specified.
  - cmap (str or plt.cm): which colormap to use. Default: 'Reds'
  - colorbar (bool): whether to add colorbar for the heatmap. Default: True
  - colorbar_orientation (str): ONLY FOR TRIANGLES, "horizontal" or "vertical". Default: "vertical"
  - x_ticks (list): a list of strings. Will be added at the bottom. THE FIRST TICK WILL BE AT THE START OF THE SIGNAL, THE LAST TICK WILL BE AT THE END.
  - fontsize (int): font size. Default: 24


# Structure Calling
 **Find A/B Compartments**
 ```config
 >>> from pyHiC.structures import AB_compartment
 >>> ab = AB_compartment(mat, n_th_eigenvector=1)
 ```
 Find A/B compartments with the input Hi-C contact map.
 Return a 1-D vector which has the same length with input map,
 sign (+ / -) indicates A or B compartment.
 - mat: (numpy.array, scipy.sparse.csr_matrix)
 - n_th_eigenvector (int): 1 or 2. Usually the 1-st eigenvector corresponds to 
 A/B compartments, but there might be some exceptions when it corresponds to two arms
 of a chromosome. If that happens, try to set this arg as 2. Default: 1
 
 **What other?**
 - TAD?
 - Loop? (High computational burden...)
 - 


# Comparing Contact Maps
 **HiCRep**
 
 ..to be done...

