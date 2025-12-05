/* Author: Adam Jones
function declarations for numerical methods library 
*/

/*
@description:
	this function implements a numerical derivative using the rk2 
numerical derivative.
	equation is of form:
		d(x(t))/dt = f(x(t), t)

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
float *euler_method(float (*x_prime)(float, float *, uint), float x[], float t[], float dt, uint N_t, uint N_x);

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
int write_to_bin(float * data, uint N, uint dims, char fname[]);
