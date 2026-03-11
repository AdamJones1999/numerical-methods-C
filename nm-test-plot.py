# C num

import numpy as np
import sys
from matplotlib import pyplot as plt
from collections.abc import Iterator, Callable

def shm_plot(data_2d: np.ndarray) -> None:
	t = data_2d[0] # indep var 
	dt = t[1]-t[0]
	r = data_2d[1] # dep var

	# analytical solution to y'' + 16y = 0
	r_analytical = np.cos(4*t) + 0.25*np.sin(4*t)

	print(f"first few: {r[0]}, {r[1]}, {r[2]}, {r[3]}, {r[4]}, {r[5]}")

	# plotting
	plt.subplot(2, 1, 1)
	plt.plot(t, r)
	plt.title(f"shmODE r over t=[{t[0]}, {t[cols-1]}], dt={dt}")
	plt.subplot(2, 1, 2)

	plt.plot(t, r-r_analytical)
	plt.title(f"shmODE r-r_analytical over t=[{t[0]}, {t[cols-1]}], dt={dt}")
	plt.show()

def holmann_orbitalmotion_test(data_2d: np.ndarray) -> None:
	nrows = data_2d.shape[0]
	ncols = data_2d.shape[1]
	if (nrows != 3):
		raise Exception(f"holmann transfer test: should have 3 rows, got {nrows}.\n")
	t = data_2d[0] # indep var 
	dt = t[1]-t[0]
	# dep vars
	x = data_2d[1]
	y = data_2d[2]
	# equation of circle to plot earth in 2d
	R_e = 6378
	theta = np.linspace(0, 2*np.pi, 100)
	earth_x = R_e * np.cos(theta)
	earth_y = R_e * np.sin(theta)

	# finding apogee
	# while orbital motion is above x axis 2 steps after start of orbit (first 180 degrees of orbit)
	i=2
	while y[i] >= 0 and i < (ncols - 1):
		i = i+1
	apogee = -x[i] - R_e
	print(f"i={i}\n")
	print(f"apogee altitude = R_e - x[{i}] = {-x[i]} - {R_e} = {apogee}\n")

	# plotting
	ax1 = plt.subplot(1, 1, 1)

	ax1.plot(earth_x, earth_y)
	ax1.plot(x, y)
	ax1.set_xlim(left=-2*R_e, right=2*R_e)
	ax1.set_ylim(bottom=-2*R_e, top=2*R_e)
	plt.title(f"orbital motion over t=[{t[0]}, {t[cols-1]}], dt={dt}")

	plt.show()

def holmann_orbitalburn_test(data_2d: np.ndarray) -> None:
	nrows = data_2d.shape[0]
	ncols = data_2d.shape[1]
	if (nrows != 3):
		raise Exception(f"holmann transfer test: should have 3 rows, got {nrows}.\n")
	t = data_2d[0] # indep var 
	dt = t[1]-t[0]
	# dep vars
	x = data_2d[1]
	y = data_2d[2]

	# equation of circle to plot earth in 2d
	R_e = 6378
	theta = np.linspace(0, 2*np.pi, 100)
	earth_x = R_e * np.cos(theta)
	earth_y = R_e * np.sin(theta)

	# finding apogee
	# while orbital motion is above x axis 2 steps after start of orbit (first 180 degrees of orbit)
	i=2
	while x[i] >= 6857.9856 and i < (ncols - 1):
		i = i+1
	print(f"i={i}\n")
	print(f"{x[i-5:i+5]}\n")

	# plotting
	ax1 = plt.subplot(1, 1, 1)

	ax1.plot(earth_x, earth_y)
	ax1.plot(x, y)
	ax1.set_xlim(left=-25000, right=25000)
	ax1.set_ylim(bottom=-25000, top=25000)
	plt.title(f"orbital burn then transfer orbit motion over t=[{t[0]}, {t[cols-1]}], dt={dt}")

	plt.show()

def run_test(test: Callable[[np.ndarray], None], data: np.ndarray):
	print(f"--------\nSTARTING plotting of: {test.__name__}\n")
	test(data)
	print(f"--------\nENDING plotting of: {test.__name__}\n")
	return;

if __name__ == "__main__":
	# read in data and dimensions
	fp = sys.argv[1] #
	rows = int(sys.argv[2]) # rows in data file
	cols = int(sys.argv[3]) # columns in data file
	data_1d = np.fromfile(fp, dtype=np.float64)

	# organize data into np arrays
	data_2d = data_1d.reshape((rows, cols))

	# @@@@@@@@@@@@ calling tests @@@@@@@@@@@@

	# run_test(shm_plot, data_2d, "shm_plot")
	# run_test(holmann_orbitalmotion_test, data_2d)
	run_test(holmann_orbitalburn_test, data_2d)
	# @@@@@@@@@@@ end of running tests @@@@@@@@@@@