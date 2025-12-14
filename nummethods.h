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
	TODO: UPDATE TO REFLECT DRDT_ARGS CHANGE
	drdt: pointer to function encoding differential equation. 
	r[]: pointer to array that solution will be put in.
	t[]: pointer to independent variable step array (t as it is usually time)
	N: number of elements in t[].
	dt: different between t and next independent var value.
	N_r: number of dependent variables in the ODE system.
@precond:
	t[]: must contain N_t elements of elements all with the same difference 
	bewteen adjacent elements.
	r[]: must contain initial value x[0] for the euler method to start with. 
		each variable in x[] must have N_t elements.
@return: 
	r: pointer to array containing solution.
*/
double *euler_method(double (*drdt)(double, double *, double *, uint), \
	double t[], double r[], double drdt_args[], double dt, uint N_t, uint N_r);

/* 
@description: 
	computes a single iteration of the euler method in place within r[], where 
	param drdt() is the derivative at t. Suitable in systems of coupled ODEs.
@params:
	drdt: func pointer to derivative function.
	t: independent variable value.
	r[]: array storing dependent variables.
	drdt_args[]: array of arguments passed to drdt as 3rd argument.
	dt: different between t and next independent var value.
	N_r: number of dependent variables in the ODE system.
@precond:
	r[]: must not point to last element in array or last element in row when 
	there are >1 dependent variables or undefined behaviour will happen due 
	to this euler method being a type of forward difference method.
@return:
	r: pointer to array containing solution.
*/
double *euler_single(double (*drdt)(double, double *, double *, uint), \
	double t, double r[], double drdt_args[], double dt, uint N_r);

//currently only single var
void rk4_single(double (*drdt)(double, double *, double *, uint), \
	double t, double r[], double drdt_args[], double dt, uint N_r);

//currently only single var tested
double *midpoint_single(double (*drdt)(double, double *, double *, uint), \
	double t, double r[], double drdt_args[], double dt, uint N_r);


/*
@description:
	currently only single var!!
	writes 1D array data (2D support coming soon) to a 
little endian binary file.

@params:
	data: pointer to data array.
	dims: number of dimensions in data.
	fname: string of desired name of file.

@return:
	0 if successful, -1 if not. 
*/
int to_bin(double * data, uint N, uint dims, char fname[], char mode[]);
