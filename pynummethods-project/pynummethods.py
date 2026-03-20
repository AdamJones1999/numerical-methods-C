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
def jacobian(f: Callable[[np.ndarray, np.ndarray, int], None], jacob: np.ndarray, r: np.ndarray, delta: float, Nr: int) -> None:
	f_r = np.ndarray(Nr)
	f(f_r, r, Nr) # compute <f(<r>)>
	f_perturbed = np.ndarray(Nr)
	r_perturbed = r.copy() # shallow copy is enough for floats (immutable obj)
	# compute one column of jacobian per loop iter
	for i in range(0, Nr):
		r_perturbed[i] += delta
		f(f_perturbed, r_perturbed, Nr)
		jacob[:, i] = (f_perturbed - f_r) / delta
		r_perturbed[i] -= delta # unperturb so only one var is perturbed at a time.

