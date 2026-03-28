"""
Collection of functions representing linear systems of equationse system. 
They are to be used within certain numerical methods defined in pynummethods.py 
such as the jacobian
"""

import numpy as np

"""ODE system for testing newton rhapson multivar method"""
def NR_test_ls(f_r: np.ndarray, r: np.ndarray, Nr: int):
	nvars = 2
	if (not (f_r.shape[0] == r.shape[0] == Nr) and Nr != 2):
		print(f"ERROR: # of elements in f_r and r must equal Nr and Nr must = {nvars}")
	x1 = r[0]
	x2 = r[1]
	f_r[0] = x2 * (1 + x1**2) - np.exp(-np.exp(x1+x2))
	f_r[1] = -x2 * np.sin(x1) + 0.5 - x1 * np.cos(x2)
	return