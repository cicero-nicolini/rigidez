import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection

lines = [
	[(0, 1), (1, 1)],
	[(2, 3), (3, 3)],
	[(1, 2), (1, 3)],
]

xpoints = np.array([1, 8])
ypoints = np.array([3, 10])


lc = LineCollection(lines, linewidths=2)
fig, ax = plt.subplots()
ax.add_collection(lc)
ax.plot(xpoints, ypoints, 'ro')
ax.autoscale()
ax.margins(0.1)
plt.show()