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

void Jac_test_ls(double f_r[], double r[], uint Nr) {
	uint nvars = 4;
	if (Nr != nvars) {
		printf("ERROR: # of elements in f_r and r must equal Nr and Nr must = %d\n", nvars);
		exit(1);
	}
	double x1 = r[Nr-4];
	double x2 = r[Nr-3];
	double x3 = r[Nr-2];
	double x4 = r[Nr-1];

	f_r[0] = 4.0*x1 +5.0*x2 -2.0*x3 +2.0*x4 - 6.0;
	f_r[1] = -1.0*x1 -4.0*x2 +6.0*x3 +3.0*x4 - 1.0;
	f_r[2] = -5.0*x1 +8.0*x2 -6.0*x3 +8.0*x4 - 5.0;
	f_r[3] =  1.0*x1 +2.0*x2 +3.0*x3 +4.0*x4 - 1.0;
	return;
}
