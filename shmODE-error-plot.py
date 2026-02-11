import numpy as np
import sys
from matplotlib import pyplot as plt

# read in data and dimensions
fp = sys.argv[1] #


#method_name = sys.argv[2]
rows = int(sys.argv[2]) # rows in data file
cols = int(sys.argv[3]) # columns data file
data_1d = np.fromfile(fp, dtype=np.float64)

# organize data into np arrays
data_2d = data_1d.reshape((rows, cols))
t = data_2d[0] # indep var 
dt = t[1]-t[0]
r = data_2d[1] # dep var

# analytical solution to y'' + 16y = 0
r_analytical = np.cos(4*t) + 0.25*np.sin(4*t)

print(f"first few: {r[0]}, {r[1]}, {r[2]}, {r[3]}, {r[4]}, {r[5]}")

# plotting
plt.subplot(2, 1, 1)
plt.plot(t, r)
plt.title(f"shmODE r over t=[{t[0]}, {t[cols-1]}], dt={dt}")
plt.subplot(2, 1, 2)

plt.plot(t, r-r_analytical)
plt.title(f"shmODE r-r_analytical over t=[{t[0]}, {t[cols-1]}], dt={dt}")
plt.show()
