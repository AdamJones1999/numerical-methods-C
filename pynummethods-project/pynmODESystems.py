"""
Collection of functions representing ODE systems that return slopes 
for the derivative of each 1st order ODE in the system. 
They are to be used within numerical methods defined in pynummethods.py
"""

import numpy as np
import pynummethods as nm
from matplotlib import pyplot  as plt

def Nvars_err_msg(f_name: str, Nvars: int, Nvars_req: int) -> str:
	return f"{f_name} requires {Nvars_req} dep vars. Provided {Nvars}"

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


"""
1D time independent shrodinger equation
x: independent variable
rcol: psi and phi (d(psi)dx)
Nr: number of vars in rcol
"""
def Schrodinger1DFixedE(x, rcol, Nr):
	Nr_req = 2
	if (Nr != Nr_req):
		raise Exception(Nvars_err_msg("Schrodinger1DFixedE", Nr, Nr_req))
	else:
		drdt = np.zeros((Nr), dtype=float)
		h_bar = 6.582e-16 # [eV*s] reduced planck constant
		m = 0.511e6 / (3e8) ** 2 # [eV/(m/s)^2]
		V = 0 # potential
		E = 17e3 # [eV] particle kinetic energy
		drdt[0] = rcol[1]
		drdt[1] = 2 * m/h_bar**2 * (V - E) * rcol[0]
		return drdt

"""
1D time independent shrodinger equation with E as an element in rcol.
ONLY to be used within Schrodinger1D_boundary_val
x: independent variable
rcol: [0]: psi, [1]: phi (d(psi)dx), [2]: Kinetic energy 
Nr: number of vars in rcol
"""
def Schrodinger1DVariableE(x, rcol, Nr):
	Nr_req = 3
	if (Nr != Nr_req):
		raise Exception(Nvars_err_msg("Schrodinger1DVariableE", Nr, Nr_req))
	else:
		drdt = np.zeros((Nr), dtype=float)
		h_bar = 6.582119569e-16 # [eV*s] reduced planck constant
		m = 0.51099895069e6 / (2.99792458e8) ** 2 # [eV/(m/s)^2]
		V = 0 # potential
		E = rcol[2]
		drdt[0] = rcol[1]
		drdt[1] = 2 * m/h_bar**2 * (V - E) * rcol[0]
		drdt[2] = 0
		return drdt

"""
r_last: the last element of the rk4 propagation
E: 1x1 ndarray containing only E: the kinetic energy guess [eV]
Nr: number of variables that should be 1 because we are concerned with just the initial energy
"""
def Schrodinger1D_boundary_val(r_last: np.ndarray, E: np.ndarray, Nr: int):
	Nr_req = 1
	if (Nr != Nr_req):
		raise Exception(Nvars_err_msg("Schrodinger1D_boundary_val", Nr, Nr_req))
	else:
		# initialization for rk4
		x0 = 0
		dx = 1e-14
		xf = 1e-11 # length of potential well L
		Nx = int(np.floor(1000*(xf-x0)/(1000*dx))) # Nx=1000 without scale up (too small numbers)
		# print(f"Nx: {Nx}")
		x = np.linspace(x0, xf, num=Nx, endpoint=False, dtype=float)
		r = np.zeros((Nr+2, Nx), dtype=float)
		r0 = [0, 1e-3, E[0]]
		r[:, 0] = r0
		# propagating solution attempt using rk4
		for i in range(0, Nx-1):
			nm.rk4_single(Schrodinger1DVariableE, x[i], r, i, dx, Nr+2)
		r_last[0] = r[0, Nx-1]
		return

