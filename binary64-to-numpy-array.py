import numpy as np
import sys
from matplotlib import pyplot as plt

#fp = "data/test1.data"
fp = sys.argv[1] #
rows = int(sys.argv[2]) # rows in data file
cols = int(sys.argv[3]) # columns data file
print(fp)

data_1d = np.fromfile(fp, dtype=np.float64)

data_2d = data_1d.reshape((rows, cols))
t = data_2d[0] # indep var 
r = data_2d[1] # dep var

print("t: ")
print(t)
print("r: ")
print(r)
# ======== plotting ========
plt.plot(t, r)
plt.show()

'''
plt.plot(r, )
ax
'''