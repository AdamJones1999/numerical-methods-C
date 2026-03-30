import numpy as np
from collections.abc import Iterator, Callable

def euler_single(drdt: Callable[[float, np.ndarray, int], np.ndarray], t: float, r: np.ndarray, i: int, dt: float, Nr: int) -> None:
	r[:, i+1] = r[:, i] + drdt(t, r[:, i], Nr) * dt
	return


def midpoint_single(drdt: Callable[[float, np.ndarray, int], np.ndarray], t: float, r: np.ndarray, i: int, dt: float, Nr: int) -> None:
	dt_mp = 0.5 * dt
	euler_single(drdt, t, r, i, dt_mp, Nr)
	r[:, i+1] = r[:, i] + drdt(t + dt_mp, r[:, i+1], Nr) * dt
	return


def rk4_single(drdt: Callable[[float, np.ndarray, int], np.ndarray], t: float, r: np.ndarray, i: int, dt: float, Nr: int) -> None:
	dt_mp = 0.5 * dt
	t_mp = t + dt_mp
	k1 = drdt(t, r[:, i], Nr) # slope at current point
	k2 = drdt(t_mp, r[:, i] + k1 * dt_mp, Nr) # slope at midpoint using k1
	k3 = drdt(t_mp, r[:, i] + k2 * dt_mp, Nr) # slope at midpoint using k2
	k4 = drdt(t + dt, r[:, i] + k3 * dt, Nr) # slope at endpoint using k3
	k = k1 / 6.0 + k2 / 3.0 + k3 / 3.0 + k4 / 6.0 # weighted avg of k1,2,3,4
	r[:, i+1] = r[:, i] + k * dt # rk4: forward difference using combined slope k
	return


# compute jacobian of vector valued multivariable function <f(<r>)>
# result is written to 'jacob' an NrxNr matrix
def jacobian(f: Callable[[np.ndarray, np.ndarray, int], None], jacob: np.ndarray, \
	r: np.ndarray, perturb: float, Nr: int) -> None:
	f_r = np.ndarray(Nr)
	f(f_r, r, Nr) # compute <f(<r>)>
	f_perturbed = np.ndarray(Nr)
	r_perturbed = r.copy() # shallow copy is enough for floats (immutable obj)
	# compute one column of jacobian per loop iter
	for i in range(0, Nr):
		r_perturbed[i] += perturb
		f(f_perturbed, r_perturbed, Nr)
		jacob[:, i] = (f_perturbed - f_r) / perturb
		r_perturbed[i] -= perturb # unperturb so only one var is perturbed at a time.

"""
@description:
	uses newton rhapson method to find the roots of a system of multivariable nonlinear equations represented by vector of functions <f>. This function finds a solution vector, not a particular solution vector.
@params:
	f: vector of functions
	r0: initial roots guess vector.
	target: amount that roots must change less than between iterations to consider the solution converged. 
	Nr: number of equations and number of variables (must be equal)
@return:
	r_guess: found solution vector.
"""
def newton_rhapson(f: Callable[[np.ndarray, np.ndarray, int], None], r0: np.ndarray, target,\
	Nr: int) -> np.ndarray:
	f_r = np.ndarray(shape=(Nr,), dtype=np.float64)
	jacob = np.ndarray((Nr, Nr), dtype=np.float64)
	r_converged = np.ndarray((Nr,), dtype=np.float64)
	r_guess = r0.copy() # initial NR guess
	perturb = target * 1e2 # jacobian dep var step
	converged = False # flag of if all root vector elements converge (change less than 'target')
	n_iters = 0 # performance tracking
	while (not converged):
		converged = True
		# compute <f(<r>)>
		f(f_r, r_guess, Nr)
		# compute jacobian
		jacobian(f, jacob, r_guess, perturb, Nr)
		# solve A*<dx> = <y> where A = jacobian of <f>, <dx> = <r_guess> - <true roots>, <y> = <f(<r_guess>)>
		dx = np.linalg.solve(jacob, f_r)
		# determine if all roots have converged
		i = 0
		while (converged and i < Nr):
			# if r_guess not converged, compute new r_guess for next iter
			if (abs(dx[i]) > target):
				converged = False
				r_guess = r_guess - dx
			i += 1
		n_iters += 1 # performance tracking
	print(f"Number of newton rhapson iternations needed for convergence: {n_iters}") # performance tracking
	return r_guess





