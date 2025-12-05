#include <stdlib.h>
#include <stdio.h>
#include <nummethods.h>

/*
@description: 
	simple dy/dt calculator for ODE y = y' for certain value t and r(t). 
Number of vars in N_r = 1 in this case but is still contained in a 
vector for compatibility with multivariable numerical methods that use pointers 
to arrays as arguments.
@params:
	t: value of independent variable t.
	r[]: vector of dependent var r_i values (should be 1 for this ODE)
	N_r: number of independent variables.
@return:
	k: calculated dr/dt.
*/
float simpleODE1(float t, float r[], uint N_r) {
	float k; // dr/dt
	(void)t; // cast to void to suppresses unused param warning
	if (N_r != 1) {
		fprintf(stderr, "ERROR: Given %d dependent var values, should be 1.", N_r);
		exit(1);
	}
	else {
		k = r[0];
		return k;
	}
}

int main() {
	// ======== testing euler_method() ========
	fprintf(stdout, "IN MAIN\n\n");
	char *fn = "test.data";
	uint N_t = 5;
	uint N_r = 1;
	float *t = (float *) malloc(N_t * sizeof(float));
	float *r = (float *) malloc(N_t * N_r * sizeof(float));
	uint i;
	for (i=0; i<N_t; i++) {
		t[i] = i+1;
	}
	r[0] = 1;
	float dt = 1;
	euler_method(simpleODE1, r, t, dt, N_t, N_r);
	printf("%d element solution r: \n{ %f, %f, %f, %f, %f }\n", N_t, r[0], r[1], r[2], r[3], r[4]);
	if (write_to_bin(r, 5, 1, fn) == 0) {
		free(t);
		free(r);
	}
	else {
		printf("writing to %s failed\n", fn);
	return 0;
	}
}
