#include <stdlib.h>
#include <stdio.h>

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
	float dummy = t; // suppresses unused param warning
	if (N_r != 1) {
		fprintf( stderr, "ERROR: Given %d dependent var values, should be 1.", N_r);
		exit(1);
	}
	else {
		k = r[0];
		return k;
	}
}

int main() {
	return 0;
}