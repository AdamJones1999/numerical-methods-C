#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <nummethods.h>
#include <cblas.h>
#include <lapacke.h>
#include <pthread.h>

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
void euler_helper(double r[], double r_next[], double drdt[], \
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
	double dt_mp = 0.5 * dt;
	double t_mp = t + dt_mp;
	double *r_curr = (double *) malloc(Nr * sizeof(double));
	double **k = alloc_2d_array(5, Nr); // array of intermediate rk4 slopes
	double **r_k = alloc_2d_array(3, Nr); // array of intermediate r propagations required for rk4
	double *k1 = k[0]; // first of four slopes used in weighted rk4 slope k
	double *k2 = k[1];
	double *k3 = k[2];
	double *k4 = k[3];
	double *k_weighted = k[4]; // weighted sum of k1,2,3,4
	double *r_mp_k1 = r_k[0]; // propagation of r by 0.5*dt using k1 as slope
	double *r_mp_k2 = r_k[1]; // propagation of r by 0.5*dt using k2 as slope
	double *r_ep_k3 = r_k[2]; // propagation of r by dt using k3 as slope
	uint i;
	// read r column into r_curr array
	for (i = 0; i < Nr; i++) {
		r_curr[i] = r[i][j_r];
	}
	// calculate k1 slopes
	drdt_f(k1, r_curr, t, Nr);
	// propagate r to midpoint using k1 slopes
	euler_helper(r_curr, r_mp_k1, k1, dt_mp, Nr);
	// calculate k2 slopes at midpoint
	drdt_f(k2, r_mp_k1, t_mp, Nr);
	// propagate r to midpoint using k2 slopes
	euler_helper(r_curr, r_mp_k2, k2, dt_mp, Nr);
	// calculate k3 slopes at midpoint
	drdt_f(k3, r_mp_k2, t_mp, Nr);
	// propagate r to endpoint using k3 slopes
	euler_helper(r_curr, r_ep_k3, k3, dt, Nr);
	// calculate k4 slopes at endpoint
	drdt_f(k4, r_ep_k3, t + dt, Nr);
	// calculate the weighted slope then propagate r to next solution step
	for (i = 0; i < Nr; i++) {
		k_weighted[i] = k1[i] / 6.0 + k2[i] / 3.0 + k3[i] / 3.0 + k4[i] / 6.0;
		r[i][j_r + 1] = r[i][j_r] + k_weighted[i] * dt;
	}
	free(r_curr);
	free_2d_array((void **) k);
	free_2d_array((void **) r_k);
	return r;
}


double **rk4_single_2(void (*drdt_f)(double *, double *, double, uint, void *), \
	double t, double **r, uint j_r, double dt, uint Nr, void *params) {
	double dt_mp = 0.5 * dt;
	double t_mp = t + dt_mp;
	double *r_curr = (double *) malloc(Nr * sizeof(double));
	double **k = alloc_2d_array(5, Nr); // array of intermediate rk4 slopes
	double **r_k = alloc_2d_array(3, Nr); // array of intermediate r propagations required for rk4
	double *k1 = k[0]; // first of four slopes used in weighted rk4 slope k
	double *k2 = k[1];
	double *k3 = k[2];
	double *k4 = k[3];
	double *k_weighted = k[4]; // weighted sum of k1,2,3,4
	double *r_mp_k1 = r_k[0]; // propagation of r by 0.5*dt using k1 as slope
	double *r_mp_k2 = r_k[1]; // propagation of r by 0.5*dt using k2 as slope
	double *r_ep_k3 = r_k[2]; // propagation of r by dt using k3 as slope
	uint i;
	// read r column into r_curr array
	for (i = 0; i < Nr; i++) {
		r_curr[i] = r[i][j_r];
	}
	// calculate k1 slopes
	drdt_f(k1, r_curr, t, Nr, params);
	// propagate r to midpoint using k1 slopes
	euler_helper(r_curr, r_mp_k1, k1, dt_mp, Nr);
	// calculate k2 slopes at midpoint
	drdt_f(k2, r_mp_k1, t_mp, Nr, params);
	// propagate r to midpoint using k2 slopes
	euler_helper(r_curr, r_mp_k2, k2, dt_mp, Nr);
	// calculate k3 slopes at midpoint
	drdt_f(k3, r_mp_k2, t_mp, Nr, params);
	// propagate r to endpoint using k3 slopes
	euler_helper(r_curr, r_ep_k3, k3, dt, Nr);
	// calculate k4 slopes at endpoint
	drdt_f(k4, r_ep_k3, t + dt, Nr, params);
	// calculate the weighted slope then propagate r to next solution step
	for (i = 0; i < Nr; i++) {
		k_weighted[i] = k1[i] / 6.0 + k2[i] / 3.0 + k3[i] / 3.0 + k4[i] / 6.0;
		r[i][j_r + 1] = r[i][j_r] + k_weighted[i] * dt;
	}
	free(r_curr);
	free_2d_array((void **) k);
	free_2d_array((void **) r_k);
	return r;
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


double **jacobian(void (*f)(double *, double *, uint), double **jacob_mat, double r[], \
	double perturb, uint Nr) {
	double *f_r = malloc(Nr * sizeof(double));
	double *f_perturbed = malloc(Nr * sizeof(double));
	double *r_perturbed = malloc(Nr * sizeof(double));
	memcpy(r_perturbed, r, (size_t) Nr * sizeof(double));
	f(f_r, r, Nr);
	uint j; uint i;
	for (j = 0; j < Nr; j++) {
		r_perturbed[j] += perturb; // perturb a single variable
		f(f_perturbed, r_perturbed, Nr);
		for (i = 0; i < Nr; i++) {
			jacob_mat[i][j] = (f_perturbed[i] - f_r[i]) / perturb;
		}
		r_perturbed[j] -= perturb;
	}

	free(f_r);
	free(f_perturbed);
	free(r_perturbed);
	return jacob_mat;
}

/*
calculates a column of jacobian in place
*/
void *jacobian_mt_calc_column(void *jacobian_data) {
	Jacmt_data *jacdata_p = (Jacmt_data *) jacobian_data;
	double *r_perturbed = malloc(jacdata_p->Nr * sizeof(double));
	double *f_perturbed = malloc(jacdata_p->Nr * sizeof(double));

	r_perturbed[jacdata_p->j] += jacdata_p->perturb; // perturb a single variable
	jacdata_p->f(f_perturbed, r_perturbed, jacdata_p->Nr); // compute <f_perturbed> where one var r[j] in <r> is perturbed
	// compute each element in column of jacobian
	uint i;
	for (i = 0; i < jacdata_p->Nr; i++) {
		jacdata_p->jacob[i][jacdata_p->j] = (f_perturbed[i] - jacdata_p->f_r[i]) / jacdata_p->perturb;
	}
	free(r_perturbed);
	free(f_perturbed);
	return NULL;
}


double **jacobian_mt(void (*f)(double *, double *, uint), double **jacob_mat, double r[], \
	double perturb, uint Nr) {
	double *f_r = malloc(Nr * sizeof(double)); 
	double *f_perturbed = malloc(Nr * sizeof(double)); // <f> evaled with one el in <r> perturbed
	double *r_perturbed = malloc(Nr * sizeof(double)); // <r>  with one el perturbed
	memcpy(r_perturbed, r, (size_t) Nr * sizeof(double));
	pthread_t *thread_ids = malloc(Nr * sizeof(pthread_t));
	Jacmt_data *thread_data = malloc(Nr * sizeof(Jacmt_data));

	// each thread must be given or create it's own r_perturbed, idk which yet.

	// initialize thread data and create threads
	f(f_r, r, Nr); // compute <f_r> unperturbed evaluation of function used by all threads.
	uint j;
	int status;
	for (j = 0; j < Nr; j++) {
		thread_data[j].f = *f;
		thread_data[j].jacob = jacob_mat;
		thread_data[j].f_r = f_r;
		thread_data[j].j = j;
		thread_data[j].perturb = perturb;
		thread_data[j].Nr = Nr;
		status = pthread_create(&thread_ids[j], NULL, jacobian_mt_calc_column, (void *) &thread_data[j]);
		if (status != 0) {
			printf("ERROR in jacobian_mt(): creation of thread %d failed with status %ld.\n", j, thread_ids[j]);
		}
	}

	// uint i;
	// loop to be replaced with loop creating threads that compute one column of jacobian each
	// for (j = 0; j < Nr; j++) {
		/* START OF work done in each thread ----------------------------------
		r_perturbed[j] += perturb; // perturb a single variable
		f(f_perturbed, r_perturbed, Nr); // compute <f_perturbed> where one var r[j] in <r> is perturbed
		// compute each element in column of jacobian
		for (i = 0; i < Nr; i++) {
			jacob_mat[i][j] = (f_perturbed[i] - f_r[i]) / perturb;
		}
		*/ // END OF work done in each thread ----------------------------------
		// r_perturbed[j] -= perturb; // not necessary in multithreaded case.
	// }

	// all threads must synchronize at this point
	free(f_r);
	free(f_perturbed);
	free(r_perturbed);
	return jacob_mat;
}


double *newton_rhapson(void (*f)(double *, double *, uint), double *r0, double target, uint Nr) {
	if (target <= 0) {
		printf("ERROR in newton_rhapson() arg: 'target' must be > 0\n");
		return NULL;
	}
	double perturb = target * 1.0e2; // jacobian fwd diff step and perturbation of each variable in jacobian.
	double *f_r = (double *) malloc(Nr * sizeof(double));
	double **jacob = (double **) alloc_2d_array(Nr, Nr); // jacobian
	double *r_guess = (double *) malloc(Nr * sizeof(double));
	double *dr = (double *) malloc(Nr * sizeof(double));
	int *ipiv = (int *) malloc(Nr * sizeof(double));
	memcpy(r_guess, r0, (size_t) Nr * sizeof(double));
	// double max_err = target + 1.0; // highest error of all roots. starts higher than target so while loop is entered.
	int info; // info for result of LAPACKE_dgesv
	uint converged = 0; // convergence flag
	uint n_iters = 0; // performance tracking
	uint i;
	uint j;
	/*
	Keep solving del(<f>)*<dx> = <f(<r_guess>)> and keep iterating until every el of <dx> < target

	NOTE: LAPACKE_dgesv overwrites <b> in <A>*<x>x = <b> with solution <x>. 
	Thus the calculated <dx> difference between <r_guess> and the true <r> will be written 
	into <f_r> each iteration.
	*/
	while (!converged) { // max_err > target
		for (i = 0; i < Nr; i++) {
			printf("r_guess[%d] = %f.\n", i, r_guess[i]);
		}
		converged = 1;
		f(f_r, r_guess, Nr); // compute <f(<r>)>
		jacobian(f, jacob, r_guess, perturb, Nr); // compute jacobian J(<f(<r>)>)
		// printf("jacob[0][0]= %f\n", jacob[0][0]);
		info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, (int) Nr, 1, jacob[0], (int) Nr, ipiv, f_r, 1); // solve del(<f>)*<dx> = <f(<r_guess>)>
		if (info < 0) {
			printf("LAPACKE_dgesv failed. info=%d.\n", info);
		}
		// for clarity while debugging. eventually dr[] can be eliminated and just use f_r[]
		for (i = 0; i < Nr; i++) {
			dr[i] = f_r[i];
		}
		i = 0;
		// check if convergence is false in any of calculated roots.
		while (converged && i < Nr) {
			if (fabs(dr[i]) > target) { // if error is greater than target, calculate new <r_guess>
				converged = 0;
				for (j = 0; j < Nr; j++) {
					r_guess[j] = r_guess[j] - dr[j];
				}
			}
			i += 1;
		}
		n_iters += 1;  // performance tracking
	}
	printf("newton rhapson method converged after %d iterations.\n", n_iters); // performance tracking
	free(f_r);
	free_2d_array((void **) jacob);
	free(dr);
	free(ipiv);
	return r_guess;
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
