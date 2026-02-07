"""
Collection of functions representing ODE systems that return slopes 
for the derivative of each 1st order ODE in the system. 
They are to be used within numerical methods defined in pynummethods.py
"""

import numpy as np

"""calculate dydt vector for y AND v for coupled ode system 
(y'=v, v' = -16y)
original ODE the system is recasted from: y''+16y = 0"""
def shm(t, rcol, Nr):
	Nr_req = 2
	if (Nr != Nr_req):
		print(f"shmODEsys requires {Nr_req} dep vars. Provided {Nr}/n")
		return
	else:
		dydt = np.zeros((2), dtype=float) # derivatives column
		dydt[0] = rcol[1] # y' = v
		dydt[1]= -16 * rcol[0] # v' = -16y
		return dydt

"""deterministic chaos ODE system"""
def rossler(t, rcol, Nr):
	Nr_req = 3
	if (Nr != Nr_req):
		print(f"shmODEsys requires {Nr_req} dep vars. Provided {Nr}/n")
		return
	else:
		drdt = np.zeros((3), dtype=float)
		a, b, c = 0.2, 0.2, 5.7
		x = rcol[0]
		y = rcol[1]
		z = rcol[2]

		drdt[0] = -y - z
		drdt[1] = x + a*y
		drdt[2] = b + z*x - c*z
		return drdt
