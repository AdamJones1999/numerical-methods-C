#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include <nummethods.h>
#include <odesystems.h>

void orbitalburn(double drdt_t[], double r_t[], double t, uint Nr) {
	(void) t; // time independent
	if (Nr != 7) {
		fprintf(stderr, "ERROR: Given %d dependent var values, should be 2.\n", Nr);
		exit(1);
	}
	// constants
	double mu = 3.986e5; // [km^3/s^2] gavitational param. calc'd from Newtons grav. const. & mass of Earth
	double g = 9.81e-3; // [km/s^2]
	double Isp = 300; // [s] specific impulse of thruster
	double T = 10; // thrust [kg*km/s^2] = 1000 * [kg*m/s^2]

	// dependent variables
	double x = r_t[0];
	double y = r_t[1];
	double z = r_t[2];
	double vx = r_t[3];
	double vy = r_t[4];
	double vz = r_t[5];
	double m = r_t[6];
	double r = sqrt(x*x + y*y + z*z);
	double v = sqrt(vx*vx + vy*vy + vz*vz);

	// 7 1st order ODEs
	drdt_t[0] = vx;
	drdt_t[1] = vy;
	drdt_t[2] = vz;
	drdt_t[3] = -mu * x / pow(r, 3) + T * vx / (m*v);
	drdt_t[4] = -mu * y / pow(r, 3) + T * vy / (m*v);
	drdt_t[5] = -mu * z / pow(r, 3) + T * vz / (m*v);
	drdt_t[6] = -T / (g*Isp);
}

/*
difference between orbitalmotion and orbitalburn is they have 
thruster off and on respectively.
*/
void orbitalmotion(double drdt_t[], double r_t[], double t, uint Nr) {
	(void) t; // time independent
	if (Nr != 7) {
		fprintf(stderr, "ERROR: Given %d dependent var values, should be 2.\n", Nr);
		exit(1);
	}
	// constants
	double mu = 3.986e5; // [km^3/s^2] gavitational param. calc'd from Newtons grav. const. & mass of Earth
	double g = 9.81e-3; // [km/s^2]
	double Isp = 300; // [s] specific impulse of thruster
	double T = 0; // thrust [kN]
	// dependent variables
	double x = r_t[0];
	double y = r_t[1];
	double z = r_t[2];
	double vx = r_t[3];
	double vy = r_t[4];
	double vz = r_t[5];
	double m = r_t[6];
	double r = sqrt(x*x + y*y + z*z);
	double v = sqrt(vx*vx + vy*vy + vz*vz);

	// 7 1st order ODEs
	drdt_t[0] = vx;
	drdt_t[1] = vy;
	drdt_t[2] = vz;
	drdt_t[3] = -mu * x / pow(r, 3) + T * vx / (m*v);
	drdt_t[4] = -mu * y / pow(r, 3) + T * vy / (m*v);
	drdt_t[5] = -mu * z / pow(r, 3) + T * vz / (m*v);
	drdt_t[6] = -T / (g*Isp);
}

/*
1D time independent shrodinger equation with E as an element in rcol.
ONLY to be used within Schrodinger1D_boundary_val
x: independent variable
rcol: [0]: psi, [1]: phi (d(psi)dx), [2]: Kinetic energy 
Nr: number of vars in rcol
*/
void Schrodinger1DVariableE(double drdx[], double rcol[], double x, uint Nr, void *E) {
	(void) x; // V is in fact V(x) but V(x)=0 in whole well region 0<x<L.
	uint Nr_req = 2;
	if (Nr != Nr_req) {
		printf("ERROR in Schrodinger1DVariableE(): Nr is %d, but must be %d\n", Nr, Nr_req);
		exit(1);
	}
	double h_bar = 6.582119569e-16; // [eV*s] reduced qplanck constant
	double m = 0.51099895069e6 / 8.9875517874e16; // [m_electron/c^2] = [eV/(m/s)^2]
	double V = 0; // potential
	// double E = ((double *) E)[0]; // must know type of pointer before doing ptr arithmetic on it.
	drdx[0] = rcol[1];
	drdx[1] = 2 * m / pow(h_bar, 2) * (V - ((double *) E)[0]) * rcol[0];
	return;
}


/*
r_last: the last column of the rk4 propagation written to in place.
E: ptr to double containing only E: the kinetic energy guess [eV] for this iteration of newton rhaposn. 
Nr: number of variables that should be 1 because we are concerned with just the initial energy
*/
void Schrodinger1D_boundary_val(double r_last[], double E[], uint Nr) {
	uint Nr_req = 1;
	if (Nr != Nr_req) {
		printf("ERROR in Schrodinger1D_boundary_val(): Nr must = %d\n", Nr_req);
		exit(1);
	}
	double x0 = 0;
	double dx = 5e-17; // 5e-17 is lowest I could get with L = 1e-11 before i got a segfault after a few seconds (long time) from what I assume is a memory overrun caused by too many iterations due to some numerical precision issue involving dividing/multiplying by very small numbers.
	double xf = 1e-11;
	uint Nx = (uint) floor((xf - x0) / dx);
	double *x = (double *) malloc(Nx);
	double **r = alloc_2d_array(Nr + 1, Nx);
	uint Nr_rk4 = 2;
	// setting initial conditions
	double r0[2] = {0, 1e-3};
	r[0][0] = r0[0];
	r[1][0] = r0[1];
	// propagating solution attempt using rk4 and param for kinetic energy E.
	uint i;
	for (i = 0; i < Nx; i++) {
		rk4_single_2(Schrodinger1DVariableE, x[i], r, i, dx, Nr_rk4, (void *) E); // rk4 call here
	}
	for (i = 0; i < Nr; i++) {
		// rk4 call here
		r_last[i] = r[i][Nx - 1];
	}
	return;
} 