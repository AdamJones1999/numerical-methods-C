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
	f: pointer to function encoding differential equation. 
	x[]: pointer to array that solution will be put in.
	t[]: pointer to independent variable step array (t as it is usually time)
	N: number of elements in t[].
@precond:
	x[]: must contain initial value x[0] for the euler method to start with. 
	t[]: must contain N elements of elements all with the same difference 
	bewteen adjacent elements.
	dt: must equal the difference between each element in t.
@return: 
	x: pointer to array containing solution.
*/
float *euler_method(float (*f)(float), float x[], float t[], float dt, uint N)