#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include <nummethods.h>
#include <odesystems.h>
#include <linearsystems.h>


int malloc_2d_array_tests() {
	int pass = 0;
	double el;
	uint nrow = 5;
	uint ncol = 2;
	double **a = alloc_2d_array(nrow, ncol);
	uint i;
	uint j;
	for (i=0; i<nrow; i++) {
		for (j=0; j<ncol; j++) {
			el = (double) (i+1)*(j+1);
			a[i][j] = el;
			if (a[i][j] != el) {
				pass = -1;
			}
		}
	}
	free_2d_array((void  **) a);
	
	return pass;
} 

/*
not a real test. I don't know how to test whether accessing a data has been 
actually freed or not without causing a segfault. There is a commented line in 
this 'test' that will cause a segfault thus proving the memory has been freed
*/
int free_2d_array_tests() {
	uint nrow = 5;
	uint ncol = 2;
	double **arr = alloc_2d_array(nrow, ncol);
	print_2d_array(arr, nrow, ncol);
	free_2d_array((void  **) arr);
	printf("array should be freed if free_2d_array() is correct.\n");

	// *a should be nil but still readable as 'a' is in stack memory
	printf("Dereferencing ptr a after freeing: %p\n", (void*)*arr);
	printf("Dereferenced ptr a after freeing\n");
	// *a[0] seg faults as it is a value in heap memory trying to be read
	//printf("Dereferencing a[0]: %f\n", *a[0]);
	//printf("Dereferenced ptr a[0] after freeing\n");
	return 0;
}

/*
@description: 
	simple dy/dt calculator for ODE y = y' for certain value t and r(t). 
Number of vars in Nr = 1 in this case but is still contained in a 
vector for compatibility with multivariable numerical methods that use pointers 
to arrays as arguments.
@params:
	t: value of independent variable t.
	r[]: vector of dependent var r_i values (should be 1 for this ODE)
	Nr: number of independent variables.
@return:
	k: calculated dr/dt.
*/
double simpleODE1(double t, double r[], double drdt_args[], uint Nr) {
	(void) t;
	(void) drdt_args;
	double k; // dr/dt
	(void) t; // cast to void to suppresses unused param warning
	if (Nr != 2) {
		fprintf(stderr, "simpleODE1 error: Given %d dependent var values, should be 2.\n", Nr);
		exit(1);
	}
	else {
		k = r[0];
		return k;
	}
}

void shmODE_drdt(double drdt_t[], double r_t[], double t, uint Nr) {
	(void) t;
	if (Nr != 2) {
		fprintf(stderr, "ERROR: Given %d dependent var values, should be 2.\n", Nr);
		exit(1);
	}
	drdt_t[0] = r_t[1];
	drdt_t[1] = -16*r_t[0];
	return;
}

// //////// test 1: simpleODE1 ////////  
int euler_basic_ODE() {
	char *fn = "data/test1.data";
	uint N_t = 5;
	uint Nr = 2; // num dependent variables
	uint NDIMS = 1;
	double *t = (double *) malloc(N_t * sizeof(double)); // indep var
	t[0] = 0;
	uint i;
	double dt = 1; // timestep
	for (i=1; i<N_t; i++) { // init indep var array
		t[i] = t[i-1]+dt;
	}
	// array decl
	
	double *r = (double *) malloc(N_t * Nr * sizeof(double)); // dep vars
	// dep+indep vars output array
	// double *tr = (double *) malloc(N_t * (1 + Nr) * sizeof(double));	
	r[0] = 1; // initial dep var value for ivp
	/*
	euler_method(simpleODE1, t, r, NULL, dt, N_t, Nr);
	//printf("%d element solution r: \n{ %f, %f, %f, %f, %f }\n", N_t, r[0], r[1], r[2], r[3], r[4]);
	*/
	if (to_bin(t, N_t, NDIMS, fn, "w") == 0 && \
		to_bin(r, Nr * N_t, NDIMS, fn, "a") == 0) {
		free(t);
		free(r);
		return 0;
	}
	else {
		printf("writing to %s failed\n", fn);
		return -1;
	}
}

// //////// test: simple harmonic motion using Euler Method (smhODE) ////////  
int euler_shm() {
	char *fn = "data/test2.data";
	// initialize vars
	uint NDIMS = 1;
	uint Nr = 2;
	uint N_t = 10000;
	double *t = (double *) malloc(N_t * sizeof(double));
	t[0] = 0;
	double dt = 0.001;
	uint i;
	for (i=1; i<N_t; i++) { // init indep var array
		t[i] = t[i-1]+dt;
	}
	
	double **r = alloc_2d_array(Nr, N_t);

	// initial conditions
	r[0][0] = 1;
	r[1][0] = 1;

	uint j;
	// numerical solving loop
	for (j=0; j<N_t-1; j++) {
		euler_single(shmODE_drdt, t[i], r, j, dt, Nr);
	}

	// making 1D array for data output
	double *y = (double *) malloc(N_t * sizeof(double));
	y = r[0];
	// write output to file
	if (to_bin(t, N_t, NDIMS, fn, "w") == 0 && \
		to_bin(y, N_t, NDIMS, fn, "a") == 0) {
		free_2d_array((void **) t);
		free(r);
		return 0;
	}
	else { 
		printf("writing to %s failed\n", fn);
		return -1;
	}
}

// //////// test: simple harmonic motion using Midpoint Method (smhODE) ////////  
int midpoint_shm() {
	char *fn = "data/test3.data";
	uint NDIMS = 1;
	uint N_t = 10000;
	double dt = 0.001;
	uint Nr = 2;
	double *t = (double *) malloc(N_t * sizeof(double));
	double **r = alloc_2d_array(Nr, N_t);
	// initial conditions
	r[0][0] = 1;
	r[1][0] = 1;
	t[0] = 0;
	uint i;
	for (i=1; i<N_t; i++) { // init indep var array
		t[i] = t[i-1]+dt;
	}
	// numerical solving loop
	uint j;
	for (j=0; j<N_t-1; j++) {
		//drdt_args[0] = drdt[i];
		midpoint_single(shmODE_drdt, t[j], r, j, dt, Nr);
	}

	// making 1D array for data output
	double *y = (double *) malloc(N_t * sizeof(double));
	y = *r;
	// write output to file
	if (to_bin(t, N_t, NDIMS, fn, "w") == 0 && \
		to_bin(y, N_t, NDIMS, fn, "a") == 0) {
		free(t);
		free_2d_array((void **) r);
		return 0;
	}
	else { 
		printf("writing to %s failed\n", fn);
		return -1;
	}
}

int rk4_test() {
	char *fn = "data/test4.data";
	uint NDIMS = 1;
	uint Nt = 10000;
	double dt = 0.001;
	uint Nr = 2;
	/* acceptable error slope for var in r[0] (y, in the ODE) 
	found via heuristic of running rk4 simulation
	and reading the rate of change of error off the graph and rounding up.
	units: [dependent var/indep var]*/
	double r0_err_slope = 1e-11; // 8.79 *10^(-12) actual slope
	double *t = (double *) malloc(Nt * sizeof(double));
	double **r = alloc_2d_array(Nr, Nt);
	double *r_analytical = (double *) malloc(Nt * sizeof(double));
	// initial conditions
	r[0][0] = 1;
	r[1][0] = 1;
	t[0] = 0;
	uint i;
	for (i=1; i<Nt; i++) { // init indep var array
		t[i] = t[i-1]+dt;
	}
	// error analysis
	//double err_max = 0;
	//double i_err_max = 0;

	/* loop to solve numerically and analytically and 
	verify the numerical solution is within 
	heuristically determined error bounds */
	uint j;
	for (j = 0; j < Nt; j++) {
		rk4_single(shmODE_drdt, t[j], r, j, dt, Nr);
		// analytical solution to y'' + 16y = 0
		r_analytical[j] = cos(4.0*t[j]) + 0.25 * sin(4.0*t[j]);
		if ((r[0][j] - r_analytical[j]) > (r0_err_slope * t[j])) {
			printf("ERROR: element r[0][%d] = %E has error %E > acceptable error %E.\n", \
				j, r[0][j], r[0][j] - r_analytical[j], r0_err_slope * t[j]);
			return -1;
		}
	}

	// making 1D array for data output
	double *y = (double *) malloc(Nr * Nt * sizeof(double));
	y = *r;
	// write output to file
	if (to_bin(t, Nt, NDIMS, fn, "w") == 0 && \
		to_bin(y, Nt, NDIMS, fn, "a") == 0) {
		free(t);
		free_2d_array((void **) r);
		return 0;
	}
	else { 
		printf("writing to %s failed\n", fn);
		return -1;
	}
}

/*
origin of xyz coord system is center of the earth at equator.
Earth is assumed to be a sphere.
*/
int rk4_orbitalmotion_test() {
	char *fn = "data/test5.data";
	uint NDIMS = 1;
	uint Nt = 100*60*100; // 100 minutes for dt = 0.01s
	double dt = 0.01; // [s]
	uint Nr = 7;
	double perigee = 480; // [km] perigee of Low Earth Orbit (LEO)
	double R_e = 6378; // [km] radius of earth
	double m0 = 2000; // initial mass of rocket
	double apogee; // [km]
	double target_apogee_err = 0.05; // [km] error tolerance for apogee target of 800km
	double *t = (double *) malloc(Nt * sizeof(double));
	double **r = alloc_2d_array(Nr, Nt);
	/*initial conditions
	row 0 to 2: init xyz coords, 3 to 5: init xyz velocities, 6: init mass
	*/
	r[0][0] = R_e + perigee; // [km] x position
	r[1][0] = 0; // [km] y position
	r[2][0] = 0; // [km] z position
	r[3][0] = 0; // [km/s] vx (x velocity)
	r[4][0] = 7.7102; // [km/s] vy
	r[5][0] = 0; // [km/s] vz
	r[6][0] = m0; // [kg] mass of rocket
	t[0] = 0;
	uint i;
	for (i=1; i < Nt; i++) { // init indep var array
		t[i] = t[i-1]+dt;
	}

	uint j;
	for (j = 0; j < Nt - 1; j++) {
		rk4_single(orbitalmotion, t[j], r, j, dt, Nr);
	}

	// verifying apogee of 800 km
	j = 2;
	// start 2 steps after orbit begins and trigger when orbit crosses x axis (y = 0)
	while (r[1][j] >= 0.0 && j < Nt - 1) {
		j++;
	}
	if (j < Nt - 1) {
		apogee = -r[0][j] - R_e;
		if ((apogee - 800.0) > target_apogee_err) {
			printf("ERROR: orbitalmotion() produced an apogee of %f which is more than allowed %f km away from an 800 km target.\n", apogee, target_apogee_err);
			return -1;
		}
		else {
			printf("apogee altitude found at element x[%d] = %f which is within %f of target 800 km\n", j, apogee, target_apogee_err);
		}
	}
	else {
		printf("ERROR: apogee never found\n");
		return -1;
	}

	// making 1D array for data output
	double *x = (double *) malloc(Nr * Nt * sizeof(double));
	double *y = (double *) malloc(Nr * Nt * sizeof(double));
	x = r[0];
	y = r[1];

	// write output to file
	if (to_bin(t, Nt, NDIMS, fn, "w") == 0 && \
		to_bin(x, Nt, NDIMS, fn, "a") == 0 && \
		to_bin(y, Nt, NDIMS, fn, "a") == 0) {
		free(t);
		free_2d_array((void **) r);
		return 0;
	}
	else { 
		printf("writing to %s failed\n", fn);
		return -1;
	}

}


/*
finds and displays to stdout the answers to EP428_lab2_q4d about orbital burn velocity magnitude changes and propellant mass used
*/
void ep428_lab2_q4d(double vx_end_burn, double vy_end_burn, double vz_end_burn, double v0, double propellent_used) {
	// all speed and velocities in [km/s] all masses in [kg]
	printf("@@@@@@@@@@@@@ START EP428 lab2 question 4d @@@@@@@@@@@@@\n");
	double v_end_burn = sqrt(vx_end_burn*vx_end_burn + vy_end_burn*vy_end_burn + vz_end_burn*vz_end_burn);
	double v_orbital_burn_diff = v_end_burn - v0;
	double target_v_orbital_burn_diff = 1.723; // [km/s]
	double v_orbital_burn_diff_err = v_orbital_burn_diff - target_v_orbital_burn_diff;
	printf("Rocket velocities at end of orbital burn (all in km/s):\n vx: %f\n vy: %f\n vz:  %f\n v: %f\nEnd of orbital burn speed - starting speed = %f km/s\nThis is %f km/s off of the change in velocity magnitude imparted by the thruster of %f km/s predicted by the approximate impulse theory.\n\n", \
		vx_end_burn, vy_end_burn, vz_end_burn, v_end_burn, v_orbital_burn_diff, v_orbital_burn_diff_err, target_v_orbital_burn_diff);
	printf("Amount of propellant used in orbital burn: %f kg\n", propellent_used);
	printf("@@@@@@@@@@@@@ END EP428 lab2 question 4d   @@@@@@@@@@@@@\n");
}


int rk4_orbitalburn_test() {
	char *fn = "data/test6.data";
	uint NDIMS = 1;
	double dt = 0.01; // [s]
	double total_apogee; // [km]
	double apogee_err; // [km]
	double target_total_apogee = 22378; // [km]
	// double target_apogee_err = 10; // [km] error tolerance for total apogee target
	double burn_time = 261.11; // [s]
	uint burn_end_i = (uint) round(burn_time / dt) - 1; // index of end of orbital burn time.
	uint dt_per_sec = ceil(1 / dt);
	printf("dt per sec: %d\n", dt_per_sec);
	uint Nt = burn_end_i + 250*60*dt_per_sec; // burn time + 250 minutes in Transfer Orbit
	uint Nr = 7;
	double perigee = 480; // [km] perigee of Low Earth Orbit (LEO)
	double R_e = 6378; // [km] radius of earth
	double m0 = 2000; // initial mass of rocket
	
	double *t = (double *) malloc(Nt * sizeof(double));
	double **r = alloc_2d_array(Nr, Nt);
	/*initial conditions
	row 0 to 2: init xyz coords, 3 to 5: init xyz velocities, 6: init mass
	*/
	r[0][0] = R_e + perigee; // [km] x position start at perigee of orbit
	r[1][0] = 0; // [km] y position
	r[2][0] = 0; // [km] z position
	r[3][0] = 0; // [km/s] vx (x velocity)
	r[4][0] = 7.7102; // [km/s] vy initial at perigee
	r[5][0] = 0; // [km/s] vz
	r[6][0] = m0; // [kg] mass of rocket
	t[0] = 0;
	uint i;
	for (i=1; i < Nt; i++) { // init indep var array
		t[i] = t[i-1]+dt;
	}

	uint j;
	// orbital burn from 0s to 'burn_time' seconds to get from LEO to 
	// Transfer Ellipse Orbit
	for (j = 0; j < burn_end_i; j++) {
		rk4_single(orbitalburn, t[j], r, j, dt, Nr);
	}

	// answering ep428_lab2_q4d
	double v0 = sqrt(r[3][0]*r[3][0] + r[4][0]*r[4][0] + r[5][0]*r[5][0]);
	double orbital_burn_propellent = m0 - r[6][j-1]; // [kg]
	ep428_lab2_q4d(r[3][j-1], r[4][j-1], r[5][j-1], v0, orbital_burn_propellent);
	
	printf("index to stop orbital burn: %d\n", burn_end_i);
	printf("j after orbital burn: %d\n", j);
	// Transfer Ellipse Orbit where spacecraft can coast with no thrust, 
	// thus using orbitalmotion ODE system due to no thrust component required
	for (;j < Nt - 1; j++) {
		rk4_single(orbitalmotion, t[j], r, j, dt, Nr);
	}

	// verifying apogee of 800 km
	j = burn_end_i;
	// continue transfer orbit with booster off and trigger when orbit changes x direction
	// orbit does not necessarily turn around when y = 0. I think this is due to 
	// numerical error propagation.
	while (r[0][j] - r[0][j-1] <= 0.0 && j < Nt - 1) {
		j++;
	}
	if (j < Nt - 1) {
		total_apogee = -r[0][j];
		apogee_err = total_apogee - target_total_apogee;
		printf("found at element x[%d], orbitalburn() into orbitalmotion() produced a total apogee of \n%f which is %f km away from an 22,378 km target.\n", j, total_apogee, apogee_err);
		/*
		else {
			printf("apogee altitude found at element x[%d] = %f which is within %f of target 800 km\n", j, total_apogee, apogee_err);
		}*/
	}
	else {
		printf("ERROR: apogee never found\n");
		return -1;
	}

	// making 1D array for data output
	double *x = (double *) malloc(Nr * Nt * sizeof(double));
	double *y = (double *) malloc(Nr * Nt * sizeof(double));
	x = r[0];
	y = r[1];

	// write output to file
	if (to_bin(t, Nt, NDIMS, fn, "w") == 0 && \
		to_bin(x, Nt, NDIMS, fn, "a") == 0 && \
		to_bin(y, Nt, NDIMS, fn, "a") == 0) {
		free(t);
		free_2d_array((void **) r);
		return 0;
	}
	else { 
		printf("writing to %s failed\n", fn);
		return -1;
	}

}


int jacobian_test() {
	uint Nr = 2;
	double **jacob = alloc_2d_array(Nr, Nr);
	double correct_jacob[2][2] = {{40.04, 26.0}, {-0.50016456, -2.84136609}};
	double r0[2] = {5.0, 4.0};
	
	printf("initial root vector <r> guess r0: [%f, %f]\n", r0[0], r0[1]);
	double perturb = 0.01; // perturbation for forward difference;
	jacobian(NR_test_ls, jacob, r0, perturb, Nr);

	// verifyinng jacobian correctness
	printf("calculated jacobian:\n[ %f, %f;\n%f, %f ]\n", jacob[0][0], jacob[0][1], jacob[1][0], jacob[1][1]);
	printf("correct jacobian:\n[ %f, %f;\n%f, %f ]\n", correct_jacob[0][0], correct_jacob[0][1], correct_jacob[1][0], correct_jacob[1][1]);

	uint i; uint j;
	for (i = 0; i < Nr; i++) {
		for (j = 0; j < Nr; j++) {
			if (fabs(jacob[i][j] - correct_jacob[i][j]) > 1e-6) {
				printf("incorrect jacobian element at el [%d][%d]: %f\n", i, j, jacob[i][j]);
				return -1;
			}
		}
	}
	printf("JACOBIAN PASSED\n");
	return 0;
}


void run_test(int (*test)(), char *test_name) {
	printf("--------\nSTARTING test: %s.\n", test_name);
	int result = test();
	if (result == -1) {
		printf("FAILED test: %s\n", test_name);
	}
	else if (result == 0) {
		printf("PASSED test: %s\n", test_name);
	}
	else {
		printf("test: %s returned undefined return code of %d \
			(not 0 or -1)", test_name, result);
	}
	return;
}

int main() {
	// //////////////// START TESING ////////////////
	
	run_test(malloc_2d_array_tests, "writing then reading from 2d array with row ptrs");
	// run_test(free_2d_array_tests, "freeing array and (if uncommented) causing segfault to verify memory is freed");
	// run_test(euler_basic_ODE, "simplest 1st order ODE: y = y'");
	run_test(euler_shm, "euler method simple harmonic motion: y'' + 16y = 0");
	run_test(midpoint_shm, "midpoint method simple harmonic motion: y'' + 16y = 0");
	run_test(rk4_test, "rk4 prelim test");
	run_test(rk4_orbitalmotion_test, "rk4 hohmann ordbitalmotion model test");
	run_test(rk4_orbitalburn_test, "rk4 hohmann ordbitalburn model test");
	run_test(jacobian_test, "jacobian test on system of two non-linear equations of 2 variables");
	// //////////////// END TESING ////////////////
	return 0;
}
