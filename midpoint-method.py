import numpy as np
from numpy._core.numeric import size

"""calculate dydt vector for y AND v for coupled ode system 
(y'=v, v' = -16y)
original ODE the system is recasted from: y''+16y = 0"""
def shmODEsys(t, rcol, dt, N_r):
	dydt = np.zeros((2), dtype=float) # derivatives column
	dydt[0] = rcol[1] # y' = v
	dydt[1]= -16 * rcol[0] # v' = -16y
	print(f"rcol: {rcol}")
	print(f"rcol shape: {rcol.shape}")
	print(f"dydt: {dydt}")
	print(f"dydt shape: {dydt.shape}")
	return dydt


def euler_single(drdt, t, r, i, dt, N_r):
	print(f"r[:,{i}]: {r[:,i]}")
	print(f"r[:,{i}] shape: {r[:,i].shape}")
	r[:, i+1] = r[:, i] + drdt(t, r[:, i], dt, N_r) * dt
	return

def midpoint_single(drdt, t, r, i, dt, N_r):
	return

	#f = np.linspace(t0, tf, (tf-t0)/dt, endpoint=False)

if __name__=="__main__":
	t0 = 0
	tf = 10
	dt = 0.01
	Nt = int(np.floor((tf-t0)/dt))
	Nr = 2 # num ODE system vars
	t = np.linspace(t0, tf, num=Nt, endpoint=False, dtype=float)
	r = np.zeros((Nr, Nt), dtype=float)
	r0 = [1, 1]
	print(f"t.shape: {t.shape}")
	print(f"r.shape: {r.shape}")
	print(f"Nt shape: {Nt}")
	r[:,0] = r0
	for i in range(0,999):
		euler_single(shmODEsys, t, r, i, dt, Nr)
	print(r[0,0:5])
