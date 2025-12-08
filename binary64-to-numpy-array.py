import numpy as np

fp = "test.data"

data_arr = np.fromfile(fp, dtype=np.float64)

print(data_arr)