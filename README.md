# Commonly-used HiC Operations for Liu Lab (Python 3) #
## -- [NOT FINISHED YET] -- ##
@Author: Fan Feng

### Required Packages
- numpy
- scipy
- matplotlib
- seaborn

### Install from GitHub ###
You can install the package with following command:
  ```console
    $ git clone https://github.com/liu-bioinfo-lab/scHiCTools.git
    $ cd scHiCTools
    $ python setup.py install
  ```


### Load Hi-C Contact Maps
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
 

### Functions
 **Normalization**
 ```config
 >>> from pyHiC.normalization import normalization
 >>> normalized_mat_1 = normalization(HiC_mat, method='log', base=10)
 >>> normalized_mat_2 = normalization(HiC_mat, method='VC_SQRT')
 ```
 Normalize Hi-C contact maps, return normalized map.
 - method: 
   - "log": take log(x + 1) for each value (make sure 0 is still 0).
   Additional argument: base (default: e), base of logarithm.
   - "power": take x^p for each value (usually 0 < p < 1).
   Additional argument: power (default: 0.5), the power "p".
   - "VC": each value divided by the sum of corresponding row then
   divided by the sum of corresponding column
   - "VC_SQRT": each value divided by the sqrt of the sum of corresponding row then
   divided by the sqrt of the sum of corresponding column
   - "KR": the sum of each row / column is one
   - "OE": each value divided by the average of its corresponding strata (diagonal line)
 
 **Visualize**
 ```config
 >>> from visualization import *
 >>> visualize_one_contact_map(mat, vmax=1, save_path=None)
 >>> visualize_two_contact_maps(mat1, mat2, vmax=1, save_path='compare.png')
 ```
 Visualize one Hi-C map with heatmap or compare two Hi-C maps
 (in upper / lower triangle in the same heatmap).
 - mat / mat1 / mat2: (numpy.array, scipy.sparse.coo_matrix)
 - vmax: (int or float) maximum value for heatmap.
 - save_path: (str or None) default: None. Path for saving the figure.
 If None, it will not be saved.
 
 **Calculate Insulation Scores**
 - To be finished
 
 **Find A/B Compartments**
 ```config
 >>> from AB_compartment import AB_compartment
 >>> ab = AB_compartment(mat)
 ```
 Find A/B compartments with the input Hi-C contact map.
 Return a 1-D vector which has the same length with input map,
 sign (+ / -) indicates A or B compartment.
 - mat: (numpy.array, scipy.sparse.coo_matrix)
 
 **What other?**
 - TAD?
 - Loop? (High computation burden...)
 - 




