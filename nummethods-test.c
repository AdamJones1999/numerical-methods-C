#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <nummethods.h>


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

void shoODE_drdt(double drdt_t[], double r_t[], double t, uint Nr) {
	(void) t;
	if (Nr != 2) {
		fprintf(stderr, "ERROR: Given %d dependent var values, should be 2.\n", Nr);
		exit(1);
	}
	drdt_t[0] = r_t[1];
	drdt_t[1] = -16*r_t[0];
	return;
}

// simple harmonic motion y`` + 16y = 0 slope calc for y'
/* TODO: likely replaced by multivar shoODE_drdt()
double shoODE_dydt(double drdt_t[], double y[], double t, uint Nr) {
	(void) drdt_t;
	(void) t;
	(void) y;
	double dydt; // 1st deriv of y
	dydt = 2; // placeholder
	if (Nr != 2) {
		fprintf(stderr, "ERROR: Given %d dependent var values, should be 2.\n", Nr);
		exit(1);
	}
	else {
		// dydt = dydt_args[0];
		return dydt;
	}
}
*/
// simple harmonic motion y`` + 16y = 0 slope calc for y''
/* TODO: likely replaced by multivar shoODE_drdt()
double shoODE_dydt2(double drdt_t[], double dydt[], double t, uint Nr) {
	(void) t;
	(void) dydt;
	double dydt2; // 2nd deriv of y
	dydt2 = 34;
	if (Nr != 2) {
		fprintf(stderr, "ERROR: Given %d dependent var values, should be 2.", Nr);
		exit(1);
	}
	else {
		// dydt2 = -16*dydt2_args[0];
		return dydt2;
	}
}
*/
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
	euler_method(simpleODE1, t, r, NULL, dt, N_t, Nr);
	//printf("%d element solution r: \n{ %f, %f, %f, %f, %f }\n", N_t, r[0], r[1], r[2], r[3], r[4]);

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
int euler_sho() {
	char *fn = "data/test2.data";
	// initialize vars
	uint NDIMS = 1;
	uint Nr = 2;
	uint N_t = 1000;
	double *t = (double *) malloc(N_t * sizeof(double));
	t[0] = 0;
	double dt = 0.01;
	uint i;
	for (i=1; i<N_t; i++) { // init indep var array
		t[i] = t[i-1]+dt;
	}
	
	double **r = (double **) malloc(Nr * N_t * sizeof(double));
	// TODO: to save memory, only keep current and next dydt vals.
	//double *drdt = (double *) malloc(Nr * N_t * sizeof(double));
	// double *drdt_args = (double *) malloc(Nr * sizeof(double));
	// double *drdt2_args = (double *) malloc(Nr * sizeof(double));
	// initial conditions
	r[0][0] = 1;
	//drdt[0] = 1;
	uint j;
	// numerical solving loop
	for (j=0; j<N_t-1; j++) {
		// drdt_args[0] = drdt[i];
		euler_single(shoODE_drdt, t[i], r, j, dt, Nr);
		// drdt2_args[0] = r[i];
		//euler_single(shoODE_dydt2, t[i], &drdt[i], dt, Nr); // calc drdt[i+1]
	}

	double *y = (double *) malloc(Nr * N_t * sizeof(double));
	y = r[0];
	// write output to file
	if (to_bin(t, N_t, NDIMS, fn, "w") == 0 && \
		to_bin(y, N_t, NDIMS, fn, "a") == 0) {
		free(t);
		free(r);
		//free(drdt);
		//free(drdt_args);
		//free(drdt2_args);
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
	uint N_t = 1000;
	double dt = 0.01;
	uint Nr = 2;
	double *t = (double *) malloc(N_t * sizeof(double));
	double *r = (double *) malloc(Nr * N_t * sizeof(double));
	// TODO: to save memory, only keep current and next drdt vals.
	double *drdt = (double *) malloc(Nr * N_t * sizeof(double));
	//double *drdt_args = (double *) malloc(Nr * sizeof(double));
	//double *drdt2_args = (double *) malloc(Nr * sizeof(double));
	// initial conditions
	r[0] = 1;
	drdt[0] = 1;
	t[0] = 0;
	uint i;
	for (i=1; i<N_t; i++) { // init indep var array
		t[i] = t[i-1]+dt;
	}
	// numerical solving loop
	for (i=0; i<N_t-1; i++) {
		//drdt_args[0] = drdt[i];
		//midpoint_single(shoODE_drdt, t[i], &r[i], dt, Nr);
		//drdt2_args[0] = r[i];
		//midpoint_single(shoODE_dydt2, t[i], &drdt[i], dt, Nr); // calc drdt[i+1]
	}

	// write output to file
	if (to_bin(t, N_t, NDIMS, fn, "w") == 0 && \
		to_bin(r, N_t, NDIMS, fn, "a") == 0) {
		free(t);
		free(r);
		free(drdt);
		//free(drdt_args);
		//free(drdt2_args);
		return 0;
	}
	else { 
		printf("writing to %s failed\n", fn);
		return -1;
	}
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
	//run_test(free_2d_array_tests, "freeing array and (if uncommented) causing segfault to verify memory is freed");
	run_test(euler_basic_ODE, "simplest 1st order ODE: y = y'");
	run_test(euler_sho, "euler method simple harmonic motion: y'' + 16y = 0");
	run_test(midpoint_shm, "midpoint method simple harmonic motion: y'' + 16y = 0");
	
	// //////////////// END TESING ////////////////
	return 0;
}
