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

double smhODE(double t, double r[], uint N_rvars) {
	(void) t;
	(void) r;
	(void) N_rvars;
	double drdt = 0; // slope
	if (N_rvars != 1) {
		fprintf(stderr, "ERROR: Given %d dependent var values, should be 1.", N_rvars);
		exit(1);
	}
	else {
		;
		return drdt;
	}
}

// //////// test 1: simpleODE1 ////////  
void test1() {
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

	if (write_to_bin(t, N_t, NDIMS, fn) == 0) {
		if (write_to_bin(r, N_rvars * N_t, NDIMS, fn) == 0) {
			free(t);
			free(r);
		}
	}
	else {
		printf("writing to %s failed\n", fn);
	}
}

void test2() {
	// //////// test 2: simple harmonic motion (smhODE) ////////  
	char *fn = "data/test2.data";
	uint N_t = 1000;
	uint N_rvars = 1;
	double *t = (double *) malloc(N_t * sizeof(double));
	double *r = (double *) malloc(N_rvars * N_t * sizeof(double));
	double *drdt = (double *) malloc(N_rvars * N_t * sizeof(double));
	double dt = 0.01;
	(void) drdt;
	(void) fn;

	uint i;
	t[0] = 0;
	for (i=1; i<N_t; i++) {
		t[i] = t[i-1] + dt;
	}
	r[0] = 0;
}

int main() {
	// //////////////// testing euler_method() ////////////////

	// //////// test 1: simpleODE1 ////////  
	// var decl and init
	fprintf(stdout, "IN MAIN\n\n");
	test1();

	// //////////////// END TESING ////////////////
	return 0;
}
