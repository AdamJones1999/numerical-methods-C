#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <nummethods.h>

double *euler_method(double (*drdt)(double, double *, double *, uint), \
	double t[], double r[], double drdt_args[], double dt, uint N_t, uint N_r) {
	uint i;
	(void) drdt_args;
	// compute solution values up to index N
	for (i=0; i<N_t-1; i++) {
		r[i+1] = r[i] + drdt(t[i], &r[i], &drdt_args[i], N_r) * dt;
	} 
	return r;
}

void rk4_single(double (*drdt)(double, double *, double *, uint), \
	double t, double r[], double drdt_args[], double dt, uint N_r) {
	(void) drdt; (void) t; (void) r; (void) drdt_args; (void) dt; (void) N_r;
	return;
}

// currently only single var!!
double *euler_single(double (*drdt)(double, double *, double *, uint), \
	double t, double r[], double drdt_args[], double dt, uint N_r) {
	*(r+1) = *r + drdt(t, r, drdt_args, N_r) * dt;
	return r; 
}

double *midpoint_single(double (*drdt)(double, double *, double *, uint), \
	double t, double r[], double drdt_args[], double dt, uint N_r) {
	double dt_mp = 0.5 * dt;
	double *r_mp = euler_single( \
		drdt, t, r, drdt_args, dt_mp, N_r);
	//estimate soln using midpoint slope where midpoint r_mid.

	*(r+1) = *r + drdt(t + dt_mp, r_mp, drdt_args, N_r) * dt_mp;
	return r;
}

int to_bin(double *data, uint N, uint dims, char fname[], char mode[]) {
	FILE *fp;
	if ((fp = fopen(fname, mode)) == NULL) {
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
			printf("writing %d dimension data to binary not supported.\n", dims);
		}
		fclose(fp);
		return -1;
	}
}
