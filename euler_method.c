#include <stdlib.h>
#include <ctype.h>
#include <nummethods.h>

float *euler_method(float (*x_prime)(float, float *, uint), float x[], float t[], float dt, uint N_t, uint N_x) {
	uint i;
	// compute solution values up to index N
	for (i=0; i<N_t-1; i++) {
		x[i+1] = x[i] + x_prime(t[i], &x[i], N_x) * dt;
	} 
	return x;	
}

int to_csv(float *, dims, char name[]) {
	FILE *fp;
	if ((fp = fopen(name, 'a')) == NULL) {
		printf("cannot open file: %s", name);
	}
	else {
		if (dims == 1) {

		}
		else if (dims == 2) {
			printf("support for writing 2D data is coming soon"\n);
			return -1;
		}
		else {
			printf("writing %d dimension data to a .csv not supported.\n");
			return -1;
		}
	}
	fclose(FILE *fp);
	return 0;
}