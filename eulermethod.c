#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>

/*
this function implements a numerical derivative using the rk2 
numerical derivative.
	equation is of form:
		d(x(t))/dt = f(x(t), t)

arguments:
	f: pointer to function encoding differential equation. 
	x[]: pointer to array that solution will be put in.
	t[]: pointer to independent variable step array (t as it is usually time)
	N: number of elements in t[].
preconditions:
	x[]: must contain initial value x[0] for the euler method to start with. 
returns: 
	x: pointer to array containing solution.
*/
float *eulermethod(float (*f)(float), float x[], float t[], uint N) {

	return x;	
}
