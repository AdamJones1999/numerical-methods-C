#include <stdlib.h>
#include <stdio.h>
#include <nummethods.h>

/*
@description: 
	simple dy/dt calculator for ODE y = y' for certain value t and r(t). 
Number of vars in N_rvars = 1 in this case but is still contained in a 
vector for compatibility with multivariable numerical methods that use pointers 
to arrays as arguments.
@params:
	t: value of independent variable t.
	r[]: vector of dependent var r_i values (should be 1 for this ODE)
	N_rvars: number of independent variables.
@return:
	k: calculated dr/dt.
*/
double simpleODE1(double t, double r[], double dydt_args[], uint N_rvars) {
	(void) t;
	(void) dydt_args;
	double k; // dr/dt
	(void)t; // cast to void to suppresses unused param warning
	if (N_rvars != 1) {
		fprintf(stderr, "ERROR: Given %d dependent var values, should be 1.", N_rvars);
		exit(1);
	}
	else {
		k = r[0];
		return k;
	}
}

// simple harmonic motion y`` + 16y = 0 slope calc for y'
double shoODE_dydt(double t, double y[], double dydt_args[], uint N_rvars) {
	(void) t;
	(void) y;
	(void) dydt_args;
	double dydt ; // 1st deriv of y
	if (N_rvars != 1) {
		fprintf(stderr, "ERROR: Given %d dependent var values, should be 1.", N_rvars);
		exit(1);
	}
	else {
		dydt = dydt_args[0];
		return dydt;
	}
}

// simple harmonic motion y`` + 16y = 0 slope calc for y''
double shoODE_dydt2(double t, double dydt[], double dydt2_args[], uint N_rvars) {
	(void) t;
	(void) dydt;
	double dydt2; // 2nd deriv of y
	if (N_rvars != 1) {
		fprintf(stderr, "ERROR: Given %d dependent var values, should be 1.", N_rvars);
		exit(1);
	}
	else {
		dydt2 = -16*dydt2_args[0];
		return dydt2;
	}
}

// //////// test 1: simpleODE1 ////////  
int euler_basic_ODE() {
	char *fn = "data/test1.data";
	uint N_t = 5;
	uint N_rvars = 1; // num dependent variables
	uint NDIMS = 1;
	double dt = 1; // timestep
	// array decl
	double *t = (double *) malloc(N_t * sizeof(double)); // indep var
	double *r = (double *) malloc(N_t * N_rvars * sizeof(double)); // dep vars
	// dep+indep vars output array
	// double *tr = (double *) malloc(N_t * (1 + N_rvars) * sizeof(double));
	uint i;
	t[0] = 0;
	for (i=1; i<N_t; i++) { // init indep var array
		t[i] = t[i-1]+dt;
	}
	r[0] = 1; // initial dep var value for ivp
	euler_method(simpleODE1, t, r, NULL, dt, N_t, N_rvars);
	printf("%d element solution r: \n{ %f, %f, %f, %f, %f }\n", N_t, r[0], r[1], r[2], r[3], r[4]);

	if (write_to_bin(t, N_t, NDIMS, fn) == 0 && \
		write_to_bin(r, N_rvars * N_t, NDIMS, fn) == 0) {
		free(t);
		free(r);
		return 0;
	}
	else {
		printf("writing to %s failed\n", fn);
		return -1;
	}
}

int euler_shm() {
	// //////// test 2: simple harmonic motion (smhODE) ////////  
	char *fn = "data/test2.data";
	uint NDIMS = 1;
	uint N_t = 1000;
	uint N_yvars = 1;
	double *t = (double *) malloc(N_t * sizeof(double));
	double *y = (double *) malloc(N_yvars * N_t * sizeof(double));
	// TODO: to save memory, only keep current and next dydt vals.
	double *dydt = (double *) malloc(N_yvars * N_t * sizeof(double));
	double *dydt_args = (double *) malloc(N_yvars * sizeof(double));
	double *dydt2_args = (double *) malloc(N_yvars * sizeof(double));
	// double *dydt = (double *) malloc(N_yvars * N_t * sizeof(double));
	// initial conditions
	y[0] = 1;
	dydt[0] = 1;
	//dydt_curr_next[1] = 0;
	double dt = 0.001;
	t[0] = 0;
	uint i;
	for (i=1; i<N_t; i++) { // init indep var array
		t[i] = t[i-1]+dt;
	}
	// numerical solving loop
	for (i=0; i<N_t-1; i++) {
		dydt_args[0] = dydt[i];
		euler_single(shoODE_dydt, t[i], &y[i], dydt_args, dt, N_yvars);
		dydt2_args[0] = y[i];
		euler_single(shoODE_dydt2, t[i], &dydt[i], dydt2_args, dt, N_yvars); // calc dydt[i+1]
	}

	// write output to file
	if (write_to_bin(t, N_t, NDIMS, fn) == 0 && write_to_bin(y, N_t, NDIMS, fn) == 0) {
		free(t);
		free(y);
		free(dydt);
		free(dydt2_args);
		return 0;
	}
	else { 
		printf("writing to %s failed\n", fn);
		return -1;
	}
}

void run_test(int (*test)(), char *test_name) {
	if (test() == -1) {
		printf("test: %s failed\n", test_name);
	}
	else {
		printf("test: %s passed\n", test_name);
	}
	return;
}

int main() {
	// //////////////// testing euler_method() ////////////////

	// //////// test 1: simpleODE1 ////////  
	// var decl and init

	// //////////////// START TESING ////////////////
	
	run_test(euler_basic_ODE, "simplest 1st order ODE: y = y'");
	run_test(euler_shm, "simple harmonic motion: y'' + 16y = 0");

	// //////////////// END TESING ////////////////
	return 0;
}
