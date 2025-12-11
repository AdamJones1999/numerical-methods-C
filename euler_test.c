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
double simpleODE1(double t, double r[], uint N_rvars) {
	(void) t;
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

// slope calc for y'
double shoODE_dydt(double t, double y[], uint N_rvars) {
	(void) t;
	double dydt ; // 1st deriv of y
	if (N_rvars != 1) {
		fprintf(stderr, "ERROR: Given %d dependent var values, should be 1.", N_rvars);
		exit(1);
	}
	else {
		dydt = y[0];
		return dydt;
	}
}

// slope calc for y''
double shoODE_dydt2(double t, double y[], uint N_rvars) {
	(void) t;
	double dydt2; // 2nd deriv of y
	if (N_rvars != 1) {
		fprintf(stderr, "ERROR: Given %d dependent var values, should be 1.", N_rvars);
		exit(1);
	}
	else {
		dydt2 = -16*y[0];
		return dydt2;
	}
}

// //////// test 1: simpleODE1 ////////  
int test1() {
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

	euler_method(simpleODE1, r, t, dt, N_t, N_rvars);
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

int test2() {
	// //////// test 2: simple harmonic motion (smhODE) ////////  
	char *fn = "data/test2.data";
	uint NDIMS = 1;
	uint N_t = 1000;
	uint N_yvars = 1;
	double *t = (double *) malloc(N_t * sizeof(double));
	double *y = (double *) malloc(N_yvars * N_t * sizeof(double));
	// double *dydt = (double *) malloc(N_yvars * N_t * sizeof(double));
	// initial conditions
	y[0] = 1;
	double dydt_curr = 1;
	double dydt_next = 0;
	double dt = 0.01;
	t[0] = 0;
	uint i;
	for (i=1; i<N_t; i++) { // init indep var array
		t[i] = t[i-1]+dt;
	}
	// numerical solving loop
	for (i=0; i<N_t-1; i++) {
		y[i+1] = euler_single(shoODE_dydt, y[i], t[i], dt, N_yvars);
		dydt_next = euler_single( \
			shoODE_dydt2, dydt_curr, t[i], dt, N_yvars); // calc dydt[i+1]
		dydt_curr = dydt_next;
	}

	// write output to file
	if (write_to_bin(t, N_t, NDIMS, fn) == 0 && write_to_bin(y, N_t, NDIMS, fn) == 0) {
		free(t);
		free(y);
		return 0;
	}
	else { 
		printf("writing to %s failed\n", fn);
		return -1;
	}
}

int main() {
	// //////////////// testing euler_method() ////////////////

	// //////// test 1: simpleODE1 ////////  
	// var decl and init
	fprintf(stdout, "IN MAIN\n\n");
	test1();
	test2();

	// //////////////// END TESING ////////////////
	return 0;
}
