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
	r_guess = np.array([5.0, 4.0], dtype=float) # initial guess
	f_r = np.ndarray((2,))
	f_perturbed = np.ndarray((2,))
	r_perturbed = r_guess.copy() # shallow copy is enough for floats (immutable obj)
	
	print(f"r guess: {r_guess}")
	delta = 0.01 # perturbation for forward difference
	nm.jacobian(linsystems.NR_test_ls, jacob, r_guess, delta, Nr)
	print(f"jacob: \n{jacob}")
	# explicitly calculate jacobian
	linsystems.NR_test_ls(f_r, r_guess, Nr) # calc f_1&2
	r_perturbed[0] += delta
	linsystems.NR_test_ls(f_perturbed, r_perturbed, Nr) # calc f_perturbed_0 1&2
	r_perturbed[0] -= delta
	jacob_expl[0] = (f_perturbed - f_r) / delta
	r_perturbed[1] += delta
	linsystems.NR_test_ls(f_perturbed, r_perturbed, Nr) # calc f_perturbed_1 1&2
	r_perturbed[1] -= delta
	jacob_expl[1] = (f_perturbed - f_r) / delta
	print(f"jacob_expl: \n{jacob_expl}")

	#df1dx = linsystems.NR_test_ls(f_r, )



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