#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <nummethods.h>
#include <math.h>

void NR_test_ls(double f_r[], double r[], uint Nr) {
	uint nvars = 2;
	if (Nr != nvars) {
		printf("ERROR: # of elements in f_r and r must equal Nr and Nr must = %d\n", nvars);
		exit(1);
	}
	double x1 = r[Nr-2];
	double x2 = r[Nr-1];
	f_r[0] = x2 * (1 + pow(x1, 2)) - exp(-exp(x1+x2));
	f_r[1] = -x2 * sin(x1) + 0.5 - x1 * cos(x2);
}
