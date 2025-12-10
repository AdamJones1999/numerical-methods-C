#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <nummethods.h>


double *euler_method(double (*x_prime)(double, double *, uint), double x[], double t[], double dt, uint N_t, uint N_x) {
	uint i;
	// compute solution values up to index N
	for (i=0; i<N_t-1; i++) {
		x[i+1] = x[i] + x_prime(t[i], &x[i], N_x) * dt;
	} 
	return x;
}

double *euler_single(double (*x_prime)(double, double *, uint), double r[], double t, double dt, uint N_r) {
	uint i;
	// compute solution values up to index N
	x[i+1] = x[i] + x_prime(t[i], &x[i], N_x) * dt;
	return x;
}

int write_to_bin(double *data, uint N, uint dims, char fname[]) {
	FILE *fp;
	if ((fp = fopen(fname, "a")) == NULL) {
		printf("cannot open file: %s", fname);
		return -1;
	}
	else {
		if (dims == 1) {
			fwrite(data, sizeof(double), (size_t) N, fp);
			fclose(fp);
			return 0;
		}
		else if (dims == 2) {
			printf("support for writing 2D data is coming soon.\n");
		}
		else {
			printf("writing %d dimension data to a .csv not supported.\n", dims);
		}
		fclose(fp);
		return -1;
	}
}
