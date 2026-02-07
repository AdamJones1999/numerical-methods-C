import numpy as np

def euler_single(drdt, t, r, i, dt, Nr):
	r[:, i+1] = r[:, i] + drdt(t, r[:, i], Nr) * dt
	return

def midpoint_single(drdt, t, r, i, dt, Nr):
	dt_mp = 0.5 * dt
	euler_single(drdt, t, r, i, dt_mp, Nr)
	r[:, i+1] = r[:, i] + drdt(t + dt_mp, r[:, i+1], Nr) * dt
	return