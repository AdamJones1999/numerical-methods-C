"""
Collection of functions representing ODE systems that return slopes 
for the derivative of each 1st order ODE in the system. 
They are to be used within numerical methods defined in pynummethods.py
"""

import numpy as np

def Nvars_err_msg(f_name: str, Nvars: int, Nvars_req: int) -> str:
	return f"{f_name} requires {Nvars_req} dep vars. Provided {Nvars}/n"

def non_linear_osc(t: float, rcol: np.ndarray, Nr: int) -> np.ndarray:
	Nr_req = 2
	if (Nr != Nr_req):
		raise Exception(Nvars_err_msg("non_linear_osc", Nr, Nr_req))
	else:
		dydt = np.zeros((Nr), dtype=float)
		x = rcol[0]
		v = rcol[1]
		k = 1
		dydt[0] = v
		dydt[1] = v/t -4*k*(t**2)*x
		return dydt


"""calculate dydt vector for y AND v for coupled ode system 
(y'=v, v' = -16y)
original ODE the system is recasted from: y''+16y = 0"""
def shm(t: float, rcol: np.ndarray, Nr: int) -> np.ndarray:
	Nr_req = 2
	if (Nr != Nr_req):
		raise Exception(Nvars_err_msg("shm", Nr, Nr_req))
	else:
		dydt = np.zeros((Nr), dtype=float) # derivatives column
		dydt[0] = rcol[1] # y' = v
		dydt[1] = -16 * rcol[0] # v' = -16y
		return dydt

"""deterministic chaos ODE system"""
def rossler(t: float, rcol: np.ndarray, Nr: int) -> np.ndarray:
	Nr_req = 3
	if (Nr != Nr_req):
		raise Exception(Nvars_err_msg("rossler", Nr, Nr_req))
	else:
		drdt = np.zeros((Nr), dtype=float)
		a, b, c = 0.2, 0.2, 5.7
		x = rcol[0]
		y = rcol[1]
		z = rcol[2]

		drdt[0] = -y - z
		drdt[1] = x + a*y
		drdt[2] = b + z*x - c*z
		return drdt
