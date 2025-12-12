import numpy as np
import sys
from matplotlib import pyplot as plt

# read in data and dimensions
fp = sys.argv[1] #
rows = int(sys.argv[2]) # rows in data file
cols = int(sys.argv[3]) # columns data file
data_1d = np.fromfile(fp, dtype=np.float64)

# organize data into np arrays
data_2d = data_1d.reshape((rows, cols))
t = data_2d[0] # indep var 
r = data_2d[1] # dep var

# plotting 
plt.plot(t, r)
plt.show()
