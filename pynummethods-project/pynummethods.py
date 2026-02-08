import numpy as np

def euler_single(drdt, t, r, i, dt, Nr):
	r[:, i+1] = r[:, i] + drdt(t, r[:, i], Nr) * dt
	return

def midpoint_single(drdt, t, r, i, dt, Nr):
	dt_mp = 0.5 * dt
	euler_single(drdt, t, r, i, dt_mp, Nr)
	r[:, i+1] = r[:, i] + drdt(t + dt_mp, r[:, i+1], Nr) * dt
	return

def rk4_single(drdt, t, r, i, dt, Nr):
	dt_mp = 0.5 * dt
	t_mp = t + dt_mp
	k1 = drdt(t, r[:, i], Nr) # slope at current point
	k2 = drdt(t_mp, r[:, i] + k1 * dt_mp, Nr) # slope at midpoint using k1
	k3 = drdt(t_mp, r[:, i] + k2 * dt_mp, Nr) # slope at midpoint using k2
	k4 = drdt(t + dt, r[:, i] + k3 * dt, Nr) # slope at endpoint using k3
	k = k1 / 6.0 + k2 / 3.0 + k3 / 3.0 + k4 / 6.0 # weighted avg of k1,2,3,4
	r[:, i+1] = r[:, i] + k * dt # rk4: forward difference using combined slope k
	return
