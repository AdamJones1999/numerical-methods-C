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
	printf("addr of first element of array: %p\n", (void *) arr[0]);
	printf("addr of last element of array: %p\n", (void *) &arr[nrow-1][ncol-1]);
	uint i;
	uint j;
	for (i=0; i<nrow; i++) {
		for (j=0; j<ncol; j++)
		printf("[%d][%d] value: %f\n", i, j, arr[i][j]);
		//for (j=0; j<nrow; j++) {
		//}
	}
	printf("address pointer to row pointers points to: %p\n", (void *) arr);
	printf("address the pointer that points to row pointer 0 is stored at: %p\n", (void *) &arr);
	printf("address row pointer 0 is stored at: %p\n", (void *) &arr[0]);
	printf("address row pointer 0 points to: %p\n", (void *) arr[0]);		
	printf("address of first element in contiguous block: %p\n", (void *) &arr[0][0]);
	printf("value of first element in contiguous block accessed by 2d index: %f\n", arr[0][0]);
	printf("value of first element in contiguous block accessed by dereference: %f\n", *arr[0]);
	return;
}

/* euler method for intermediate steps in other numerical methods. */
void euler_helper(double *drdt, double t, double r[], double r_next[], \
	double dt, uint Nr) {
	uint j;
	// compute solution values up to index N
	for (j = 0; j < Nr; j++) {
		r_next[j] = r[j] + drdt[j] * dt;
	} 
	return;
}

double **rk4_single(void (*drdt_f)(double *, double *, double, uint), \
	double t, double **r, uint j_r, double dt, uint Nr) {
	(void) drdt_f; (void) dt; (void) t; (void) r; (void) dt; (void) j_r;
	double **k = alloc_2d_array(Nr, 5); // array of intermediate rk4 slopes
	double **r_k = alloc_2d_array(Nr, 3); // array of intermediate r propagations required for rk4
	double *r_curr = (double *) malloc(Nr * sizeof(ptr));
	double *r_mp_k1 = r_rk4_slopes[0]// propagation of r by 0.5*dt using k1 as slope
	double *r_mp_k2 // propagation of r by 0.5*dt using k2 as slope
	double *r_ep_k3 // propagation of r by dt using k3 as slope
	uint i;
	// read r column into r_curr array
	for (i = 0; i < Nr; i++) {
		r_curr[i] = r[i][j_r];
	}
	// get k1
	drdt_f(k[0], r_curr, t, Nr);
	for 
	
	/*
	uint j;
	double *row_ptr;
	// print using 1d indexing
	printf("\n\n%d is i\n\n", i);
	i=0;
	for (i=0; i<Nr; i++) {
		row_ptr = k[i];
		for (j=0; j<Nr; j++) {
			printf("row_ptr[%d] = %f\n", j, row_ptr[j]);
			printf("   k[%d][%d] = %f\n", i, j, k[i][j]);

		}
	}*/

	//double *r_curr = (double *) malloc(Nr * sizeof(ptr));


	return k;
}

// now multivar
double **euler_single(void (*drdt_f)(double *, double *, double, uint), \
	double t, double **r, uint j_r, double dt, uint Nr) {
	double *r_curr = (double *) malloc(Nr * sizeof(double));
	double *drdt_curr = (double *) malloc(Nr * sizeof(double));
	uint i;
	// read r column into r_curr array
	for (i = 0; i < Nr; i++) {
		r_curr[i] = r[i][j_r];
	}
	// compute drdt for r colummn and t value using r_curr so drdt_f doesn't deal with indexing.
	drdt_f(drdt_curr, r_curr, t, Nr);
	// use each value in the drdt array in calculation of next solution step
	for (i = 0; i < Nr; i++) {
		r[i][j_r+1] = r[i][j_r] + drdt_curr[i] * dt;
	}
	free(r_curr);
	free(drdt_curr);
	return r; 
}

// ERROR: euler_single returns pointer at same address passed to func.
// then r+1 is modified. Then the 2nd slope calc uses r_mp which is 
// STILL pointing to r passed into midpoint_single, meaning a single euler 
// method iteration is done. No midpoint happens.
double **midpoint_single(void (*drdt_f)(double *, double *, double, uint), \
	double t, double **r, uint j_r, double dt, uint Nr) {
	// TODO: Reimplement to support multivar. Currently broken due to euler_single()
	// and changes to how drdt() functions interface with num methods.
	double *r_curr = (double *) malloc(Nr * sizeof(double));
	double *k1 = (double *) malloc(Nr * sizeof(double));
	double *r_mp = (double *) malloc(Nr * sizeof(double));
	double *k2 = (double *) malloc(Nr * sizeof(double));
	// read r column into r_curr array
	uint i;
	for (i = 0; i < Nr; i++) {
		r_curr[i] = r[i][j_r];
	}
	double dt_mp = 0.5 * dt;

	/*
	//find slope k1 to get r_mp. 
	CANT USE EULER_SINGLE, 
	just have to implement euler here for now
	*/

	// find slope k1
	drdt_f(k1, r_curr, t, Nr);
	//find r_mp using k1
	for (i = 0; i < Nr; i++) {
		r_mp[i] = r_curr[i] + k1[i] * dt_mp; 
	}
	//get slope k2 using  r_mp as input
	drdt_f(k2, r_mp, t + dt_mp, Nr);
	for (i = 0; i < Nr; i++) {
		r[i][j_r+1] = r_curr[i] + k2[i] * dt; 
	}
	free(r_curr);
	free(k1);
	free(r_mp);
	free(k2);
	return r;

	// ERROR: euler_single returns pointer at same address passed to func.
	// then r+1 is modified. Then
	/*
	double *r_mp = euler_single( \
	drdt_f, t, r, dt_mp, Nr);
	*/
	//estimate soln using midpoint slope where midpoint r_mp is r(t+dt/2) estimated using euler.
	// WRONG: OVERWRITES THE MIDPOINT
	//*(r+1) = *r + drdt_f(t + dt_mp, r_mp, Nr) * dt;
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
