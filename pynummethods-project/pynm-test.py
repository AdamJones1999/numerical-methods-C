import numpy as np
import pynmODESystems as odes
import pynmLinearSystems as linsystems
from matplotlib import pyplot  as plt
from typing import Tuple
import pynummethods as nm


def shmODEEulerTest() -> Tuple[np.ndarray, np.ndarray, float, int]:
	t0 = 0
	tf = 10
	dt = 0.001
	Nt = int(np.floor((tf-t0)/dt))
	Nr = 2 # num ODE system vars
	t = np.linspace(t0, tf, num=Nt, endpoint=False, dtype=float)
	r = np.zeros((Nr, Nt), dtype=float)
	r0 = [1, 1]
	print(f"t.shape: {t.shape}")
	print(f"r.shape: {r.shape}")
	print(f"Nt shape: {Nt}")
	r[:,0] = r0
	for i in range(0,Nt-1):
		nm.euler_single(odes.shm, t[i], r, i, dt, Nr)
	print(f"shmODEEulerTest first few: {r[0,0:5]}")
	return t, r, dt, Nt

def shmODEMidpointTest() -> Tuple[np.ndarray, np.ndarray, float, int]:
	t0 = 0
	tf = 10
	dt = 0.001
	Nt = int(np.floor((tf-t0)/dt))
	Nr = 2 # num ODE system vars
	t = np.linspace(t0, tf, num=Nt, endpoint=False, dtype=float)
	r = np.zeros((Nr, Nt), dtype=float)
	r0 = [1, 1]
	print(f"t.shape: {t.shape}")
	print(f"r.shape: {r.shape}")
	print(f"Nt shape: {Nt}")
	r[:,0] = r0
	for i in range(0,Nt-1):
		nm.midpoint_single(odes.shm, t[i], r, i, dt, Nr)
	print(f"shmODEMidpointTest first few: {r[0,0:5]}")
	return t, r, dt, Nt


def rossler_euler_test() -> Tuple[np.ndarray, np.ndarray, float, int]:
	t0 = 0
	dt = 0.1
	tf = 500+dt
	Nt = int(np.floor((tf-t0)/dt))
	Nr = 3
	t = np.linspace(t0, tf, num=Nt, endpoint=False, dtype=float)
	r = np.zeros((Nr, Nt), dtype=float)
	r0 = [0, 1, 0]
	r[:, 0] = r0
	for i in range(0, Nt-1):
		nm.euler_single(odes.rossler, t[i], r, i, dt, Nr)
	return t, r, dt, Nt

def non_linear_osc_rk4_test() -> Tuple[np.ndarray, np.ndarray, float, int]:
	t0 = 0.1
	dt = 0.01
	tf = 30+dt
	Nt = int(np.floor((tf-t0)/dt))
	Nr = 2
	t = np.linspace(t0, tf, num=Nt, endpoint=False, dtype=float)
	r = np.zeros((Nr, Nt), dtype=float)
	r0 = [0.0099983333, 0.19999]
	r[:, 0] = r0
	for i in range(0, Nt-1):
		nm.rk4_single(odes.non_linear_osc, t[i], r, i, dt, Nr)
	return t, r, dt, Nt

def jacobian_test():
	Nr = 2
	jacob = np.ndarray((Nr, Nr))
	jacob_expl = np.ndarray((Nr, Nr))
	r0 = np.array([5.0, 4.0], dtype=float)
	f_r = np.ndarray((2,))
	f_perturbed = np.ndarray((2,))
	r_perturbed = r0.copy() # shallow copy is enough for floats (immutable obj)
	
	print(f"initial root vector <r> guess: {r0}")
	delta = 0.01 # perturbation for forward difference
	nm.jacobian(linsystems.NR_test_ls, jacob, r0, delta, Nr)

	# verifyinng jacobian correctness
	print(f"calculated jacobian:\n{jacob}")
	correct_jacob = np.array([[40.04, 26.0], [-0.50016456, -2.84136609]])
	err_list = [] # changed from -1's to err location if calculation error found
	for i in range(0, Nr):
		for j in range(0, Nr):
			if (abs(jacob[i][j] - correct_jacob[i][j]) > 1e-6):
				err_list.append((i,j))
	if (len(err_list) != 0):
		print(f"ERROR: incorrect jacobian elements: {err_list}\n")


def newton_rhapson_test():
	Nr = 2
	r0 = np.array([5.0, 4.0], dtype=float)
	target = 1e-4
	f_r = np.ndarray((Nr,))
	r_solved = nm.newton_rhapson(linsystems.NR_test_ls, r0, target, Nr)
	# verify found roots are correct
	linsystems.NR_test_ls(f_r,r_solved, Nr)
	for i in range(0, Nr):
		print(f"r_solved=[{i}] = {r_solved[i]}\nf_r_solved[{i}] = {f_r[i]}")


def Schrodinger1DFixedE_test():
	x0 = 0
	dx = 1e-14
	xf = 1e-11 # length of potential well
	Nx = int(np.floor(1000*(xf-x0)/(1000*dx))) # Nx=1000 without scale up (too small numbers)
	print(f"Nx: {Nx}")
	Nr = 2
	x = np.linspace(x0, xf, num=Nx, endpoint=False, dtype=float)
	r = np.zeros((Nr, Nx), dtype=float)
	r0 = [0, 1e-3]
	r[:, 0] = r0
	for i in range(0, Nx-1):
		nm.rk4_single(odes.Schrodinger1DFixedE, x[i], r, i, dx, Nr)
	return x, r, dx, Nx

def Schrodinger1D_boundary_val_test(E0_guess, E0_correct):
	pass
	Nr = 1
	E0 = np.array([E0_guess], dtype=np.float64) # [eV] initial guess for energy level
	target = 1e-10
	E_solved = nm.newton_rhapson(odes.Schrodinger1D_boundary_val, E0, target, Nr)
	err = abs(E_solved[0] / E0_correct - 1)
	if (err > 0.002):
		print(f"ERROR: numerical solution to energy level E: {E_solved[0]} is {err} off of correct E: {E0_correct}")
	print(f"E_solved= = {E_solved[0]}\n")

def shm_euler_midpoint_plot(t_e1, r_e1, dt_e1, Nt_e1, t_m1, r_m1, dt_m1, Nt_m1):
	# ---------- analytical solution to y'' + 16y = 0 ----------
	x_shm_analytical = np.cos(4*t_e1) + 0.25*np.sin(4*t_e1)

	# euler method
	plt.figure(figsize=(16, 8))
	plt.tight_layout()
	plt.subplot(2, 2, 1)
	plt.plot(t_e1, r_e1[0, :])
	plt.title(f"euler method shmODE r over t=[{t_e1[0]}, {t_e1[Nt_e1-1]}], dt={dt_e1}")

	plt.subplot(2, 2, 2)
	plt.plot(t_e1, r_e1[0, :] - x_shm_analytical)
	plt.title(f"euler method shmODE error over t=[{t_e1[0]}, {t_e1[Nt_e1-1]}], dt={dt_e1}")
	plt.ylabel("x-x_analytical")
	plt.xlabel("t")

	# midpoint method
	plt.subplot(2, 2, 3)
	plt.plot(t_m1, r_m1[0, :])
	plt.title(f"midpoint method shmODE using r over t=[{t_m1[0]}, {t_m1[Nt_m1-1]}], dt={dt_m1}")
	
	plt.subplot(2, 2, 4)
	plt.plot(t_m1, r_m1[0, :] - x_shm_analytical)
	plt.title(f"midpoint method shmODE r-r_analytical over t=[{t_m1[0]}, {t_m1[Nt_m1-1]}], dt={dt_m1}")
	
	plt.show()


def rossler_euler_plot(t, r, dt, Nt):
	plt.figure(figsize=(16, 8))
	plt.tight_layout()
	plt.subplot(3, 1, 1)
	plt.plot(t, r[0, :])
	plt.title(f"rosslerODEsys x vs t using euler method over t=[{t[0]}, {t[Nt-1]}], dt={dt}")
	plt.xlabel("t")
	plt.ylabel("x")

	plt.subplot(3, 1, 2)
	plt.plot(t, r[1, :])
	plt.title(f"rosslerODEsys y vs t using euler method over t=[{t[0]}, {t[Nt-1]}], dt={dt}")
	plt.xlabel("t")
	plt.ylabel("y")

	plt.subplot(3, 1, 3)
	plt.plot(t, r[2, :])
	plt.title(f"rosslerODEsys z vs t using euler method over t=[{t[0]}, {t[Nt-1]}], dt={dt}")
	plt.xlabel("t")
	plt.ylabel("z")

	plt.show()

	plt.figure(figsize=(10, 10))
	plt.tight_layout()
	plt.subplot(2, 2, 1)
	plt.plot(r[0, :], r[1, :])
	plt.title(f"rosslerODEsys x,y phase space using euler method \nover t=[{t[0]}, {t[Nt-1]}], dt={dt}")
	plt.xlabel("x")
	plt.ylabel("y")

	plt.subplot(2, 2, 2)
	plt.plot(r[0, :], r[2, :])
	plt.title(f"rosslerODEsys x,z phase space using euler method \nover t=[{t[0]}, {t[Nt-1]}], dt={dt}")
	plt.xlabel("x")
	plt.ylabel("z")

	plt.subplot(2, 2, 3)
	plt.plot(r[1, :], r[2, :])
	plt.title(f"rosslerODEsys y,z phase space using euler method \nover t=[{t[0]}, {t[Nt-1]}], dt={dt}")
	plt.xlabel("y")
	plt.ylabel("z")

	plt.show()


def non_linear_osc_rk4_plot(t, r, dt, Nt):

	# analytical soln
	k = 1 # MUST MATCH NUMERICAL SOLN k
	analytical_soln = np.sin(np.sqrt(k) * t ** 2)

	plt.figure(figsize=(16, 8))
	plt.tight_layout()
	plt.subplot(2, 1, 1)
	plt.plot(t, np.log(abs(r[0, :] - analytical_soln)))
	plt.title(f"non_linear_osc using rk4: ln(x-x_analytical) over t=[{t[0]}, {t[Nt-1]}], dt={dt}")
	plt.ylabel("ln(abs(r-r_analytical))")
	plt.xlabel("t")

	plt.show()


def Schrodinger1DFixedE_plot(x, r, dx, Nx):
	plt.figure(figsize=(16, 8))
	plt.tight_layout()
	plt.subplot(1, 1, 1)
	plt.plot(x, r[0, :])
	plt.title(f"1D time indep schrod psi(x) vs x using rk4 over 0 to L x=[{x[0]}, {x[Nx-1]}], dx={dx}")
	plt.xlabel("x")
	plt.ylabel("psi(x)")

	plt.show()

if __name__=="__main__":
	'''
	# simple harmonic oscillator ODE test
	t_e1, r_e1, dt_e1, Nt_e1 = shmODEEulerTest()
	t_m1, r_m1, dt_m1, Nt_m1 = shmODEMidpointTest()
	shm_euler_midpoint_plot(t_e1, r_e1, dt_e1, Nt_e1, t_m1, r_m1, dt_m1, Nt_m1)

	# rossler ODE system test
	t_e_r1, r_e_r1, dt_e_r1, Nt_e_r1 = rossler_euler_test()
	rossler_euler_plot(t_e_r1, r_e_r1, dt_e_r1, Nt_e_r1)
	'''
	# non linear oscillator using rk4 test
	# t_rk4_1, r_rk4_1, dt_rk4_1, Nt_rk4_1 = non_linear_osc_rk4_test()
	# non_linear_osc_rk4_plot(t_rk4_1, r_rk4_1, dt_rk4_1, Nt_rk4_1)

	# jacobian test
	jacobian_test()

	# newton rhapson test
	# newton_rhapson_test()

	# Schrodinger1DFixedE test
	# x, r, dx, Nx = Schrodinger1DFixedE_test()
	# Schrodinger1DFixedE_plot(x, r, dx, Nx)
	# Schrodinger1D_boundary_val test
	Schrodinger1D_boundary_val_test(4e3, 3760.30162) # energy level n=1
	Schrodinger1D_boundary_val_test(16e3, 15041.2065) # energy level n=2
	Schrodinger1D_boundary_val_test(34e3, 33842.7146) # energy level n=3
	Schrodinger1D_boundary_val_test(60e3, 60164.8259) # energy level n=4
	print("finito")