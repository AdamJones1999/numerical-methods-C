import numpy as np
from matplotlib import pyplot  as plt
from pynummethods import euler_single, midpoint_single

"""calculate dydt vector for y AND v for coupled ode system 
(y'=v, v' = -16y)
original ODE the system is recasted from: y''+16y = 0"""
def shmODEsys(t, rcol, Nr):
	if (Nr != 2):
		print("shmODEsys requires 2 dep vars. Provided {Nr}/n")
		return
	else:
		dydt = np.zeros((2), dtype=float) # derivatives column
		dydt[0] = rcol[1] # y' = v
		dydt[1]= -16 * rcol[0] # v' = -16y
		print(f"rcol: {rcol}")
		print(f"rcol shape: {rcol.shape}")
		print(f"dydt: {dydt}")
		print(f"dydt shape: {dydt.shape}")
		return dydt

def shmODEEulerTest():
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
	return t, r, dt, Nt

def shmODEMidpointTest():
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
		midpoint_single(shmODEsys, t, r, i, dt, Nr)
	print(r[0,0:5])
	return t, r, dt, Nt

if __name__=="__main__":
	# simple harmonic oscillator ODE test
	t_e1, r_e1, dt_e1, Nt_e1 = shmODEEulerTest()
	t_m1, r_m1, dt_m1, Nt_m1 = shmODEMidpointTest()

	# ---------- analytical solution to y'' + 16y = 0 ----------
	r_shm_analytical = np.cos(4*t_e1) + 0.25*np.sin(4*t_e1)


	# ----------------------- plotting -----------------------

	# euler method
	plt.subplot(2, 2, 1)
	plt.plot(t_e1, r_e1[0, :])
	plt.title(f"euler method shmODE r over t=[{t_e1[0]}, {t_e1[Nt_e1-1]}], dt={dt_e1}")

	plt.subplot(2, 2, 2)
	plt.plot(t_e1, r_e1[0, :] - r_shm_analytical)
	plt.title(f"euler method shmODE r-r_analytical over t=[{t_e1[0]}, {t_e1[Nt_e1-1]}], dt={dt_e1}")

	# midpoint method
	plt.subplot(2, 2, 3)
	plt.plot(t_m1, r_m1[0, :])
	plt.title(f"midpoint method shmODE using r over t=[{t_m1[0]}, {t_m1[Nt_m1-1]}], dt={dt_m1}")
	
	plt.subplot(2, 2, 4)
	plt.plot(t_m1, r_m1[0, :] - r_shm_analytical)
	plt.title(f"midpoint method shmODE r-r_analytical over t=[{t_m1[0]}, {t_m1[Nt_m1-1]}], dt={dt_m1}")
	plt.show()
