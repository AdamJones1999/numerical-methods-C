import numpy as np
import sys
from matplotlib import pyplot as plt

#fp = "data/test1.data"
fp = sys.argv[1]
print(fp)

data_1d = np.fromfile(fp, dtype=np.float64)

data_2d = data_1d.reshape((2, 5))
t = data_2d[0]
r = data_2d[1]

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