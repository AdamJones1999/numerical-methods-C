#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <nummethods.h>

// allocates one block of memory for array data and then a 
// vector of row pointers.
double **alloc_2d_array(uint nrow, uint ncol) {
	double **row_ptrs = (double **) malloc(nrow * sizeof(double *));
	double *data_block = (double *) malloc(nrow * ncol * sizeof(double));
	uint i;
	for (i=0; i<nrow; i++) {
		row_ptrs[i] = &data_block[i * ncol];
	}
	return row_ptrs;
}

void free_2d_array(void **arr) {
	free(arr[0]);
	free(arr);
	return;
}

void print_2d_array(double **arr, uint nrow, uint ncol) {
	(void) ncol;
	uint i;
	uint j;
	for (i=0; i<nrow; i++) {
		for (j=0; j<ncol; j++)
		printf("address of addr [%d][%d]: %p\n", i, j, (void *) &arr[i][j]);
		//for (j=0; j<nrow; j++) {
		//}
	}
	return;
}

// depreciated. still uses drdt_args[]
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

void rk4_single(double (*drdt)(double, double *, uint), \
	double t, double r[], double dt, uint N_r) {
	(void) drdt; (void) t; (void) r; (void) dt; (void) N_r;
	return;
}

// currently only single var!!
double *euler_single(double (*drdt)(double, double *, uint), \
	double t, double r[], double dt, uint N_r) {
	*(r+1) = *r + drdt(t, r, N_r) * dt;
	return r; 
}

// ERROR: euler_single returns pointer at same address passed to func.
// then r+1 is modified. Then the 2nd slope calc uses r_mp which is 
// STILL pointing to r passed into midpoint_single, meaning a single euler 
// method iteration is done. No midpoint happens.
double *midpoint_single(double (*drdt)(double, double *, uint), \
	double t, double r[], double dt, uint N_r) {
	double dt_mp = 0.5 * dt;
	// ERROR: euler_single returns pointer at same address passed to func.
	// then r+1 is modified. Then
	double *r_mp = euler_single( \
		drdt, t, r, dt_mp, N_r);
	//estimate soln using midpoint slope where midpoint r_mp is r(t+dt/2) estimated using euler.
	// WRONG: OVERWRITES THE MIDPOINT
	*(r+1) = *r + drdt(t + dt_mp, r_mp, N_r) * dt;
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
