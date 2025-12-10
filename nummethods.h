/* Author: Adam Jones
function declarations for numerical methods library 
*/

/*
@description:
	this function implements the euler method for solving a 1st order ODE
	ODE of this form:
		d(x(t))/dt = f(x(t), t)
	does not support solving for slopes in coupled systems 
	(slopes depend on prev slopes from other equations in system of ODEs).

@params:
	x_prime: pointer to function encoding differential equation. 
	x[]: pointer to array that solution will be put in.
	t[]: pointer to independent variable step array (t as it is usually time)
	N: number of elements in t[].
	N_r: number of vars in the ODE.
@precond:
	t[]: must contain N_t elements of elements all with the same difference 
	bewteen adjacent elements.
	x[]: must contain initial value x[0] for the euler method to start with. 
		each variable in x[] must have N_t elements.
	dt: must equal the difference between each element in t.
@return: 
	x: pointer to array containing solution.
*/
double *euler_method(double (*x_prime)(double, double *, uint), double x[], double t[], double dt, uint N_t, uint N_x);

/* 
@description: 
	computes a single iteration of the euler method. Usable for 
	systems of coupled ODEs.
	
*/
double euler_single(double (*x_prime)(double, double *, uint), double r, double t, double dt, uint N_r);

/*
@description:
	writes 1D array data (2D support coming soon) to a 
little endian binary file.

@params:
	data: pointer to data array.
	dims: number of dimensions in data.
	fname: string of desired name of file.

@return:
	0 if successful, -1 if not. 
*/
int write_to_bin(double * data, uint N, uint dims, char fname[]);
