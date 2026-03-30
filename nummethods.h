/* Author: Adam Jones
function declarations for numerical methods library 
*/

/*
@description:
	this function allocates a 2D array of doubles (TODO: test if type agnostic works) as one block of memory for the array data, then another block of memory to store the pointers to each row in the main block of memory. This function handles assigning the pointers in the pointer block to the correct addresses in the main block of memory.  

@params:
	row: number of rows in array.
	col: number of columns in array
@precond:
	system is not out of memory to allocate rows to.
	
@return: 
	rows: pointer to start of first of array.
*/
double **alloc_2d_array(uint nrow, uint ncol);

void free_2d_array(void **arr);

void print_2d_array(double **arr, uint nrow, uint ncol);

/*
************* depreciated: still uses drdt_args[] *************
@description:
	this function implements one step of the euler method for 
	solving a 1st order	ODE of this form:
		d(x(t))/dt = f(x(t), t)
	where the slopes are already found. This version of the 
	euler method is intended for intermediate steps in other numerical methods.  
@params:
	drdt: 1d array of pre-calculated slopes.
	t: independent variable value
	r[]: 1d array of var values at current soln step
	r_next[]: 1d array of var values at next soln step.
	dt: step of indep variable
	Nr: number of dependent variables.
@precond:
	r[]: every element must have a double in it.
@return: 
	r_next[]: pointer to array containing propagated solution.
*/
void euler_helper(double r[], double r_next[], double drdt[], \
	double dt, uint Nr);

/* 
TODO: update depreciated docs
@description: 
	computes a single iteration of the euler method in place within r[], where 
	param drdt() is the derivative at t. Suitable in systems of coupled ODEs.
@params:
	drdt: func pointer to derivative function.
	t: independent variable value.
	r[]: array storing dependent variable data.
	dt: different between t and next independent var value.
	N_r: number of dependent variables in the ODE system.
@precond:
	r[]: must not point to last element in array or last element in row when 
	there are >1 dependent variables or undefined behaviour will happen due 
	to this euler method being a type of forward difference method.
@return:
	r: pointer to array containing solution.
*/
double **euler_single(void (*drdt_f)(double *, double *, double, uint), \
	double t, double **r, uint j_r, double dt, uint Nr);

//currently only single var
double **rk4_single(void (*drdt_f)(double *, double *, double, uint), \
	double t, double **r, uint j_r, double dt, uint Nr);

//currently only single var tested
double **midpoint_single(void (*drdt_f)(double *, double *, double, uint), \
	double t, double **r, uint j_r, double dt, uint Nr);


double **jacobian(void (*f)(double *, double *, uint), double **jacob_mat, double r[], \
	double perturb, uint Nr);

/*
@description:
	currently only single var!!
	writes 1D array data (2D support coming soon) to a little endian binary file.

@params:
	data: pointer to data array.
	dims: number of dimensions in data.
	fname: string of desired name of file.

@return:
	0 if successful, -1 if not. 
*/
int to_bin(double * data, uint N, uint dims, char fname[], char mode[]);
