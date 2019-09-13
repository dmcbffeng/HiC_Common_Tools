from scipy.sparse import coo_matrix
import numpy as np


x = np.array([[1, 0, 0],
              [2, 1, 0],
              [0, 0, 0]])
print(isinstance(x, np.ndarray))
x = coo_matrix(x)

print(x)
x[0, 2] = 5
print(x)
