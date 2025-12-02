#include <stdlib.h>
#include <ctype.h>

float *euler_method(float (*x_prime)(float, float), float x[], float t[], float dt, uint N) {
	uint i;
	// compute solution values up to index N
	for (i=0; i<N-1; i++) {
		x[i+1] = x[i] + x_prime(t[i], x[i]) * dt;
	} 
	return x;	
}