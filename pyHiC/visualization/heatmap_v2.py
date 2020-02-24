import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Arc


s = np.random.random((100,)) + 1

fig = plt.figure()
ax = fig.add_subplot(211, aspect='auto')
# plt.fill_between(np.arange(len(s)), 0, s)
a = Arc((0.5, 0.), 20, 1, theta1=180., theta2=360.)
ax.add_patch(a)
# plt.ylim([-1, 2])
plt.show()

