#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <nummethods.h>
#include <math.h>

void hohmann_transfer(double drdt_t[], double r_t[], double t, uint Nr) {
	if (Nr != 7) {
		fprintf(stderr, "ERROR: Given %d dependent var values, should be 2.\n", Nr);
		exit(1);
	}
	// constants
	mu = 3.986e8; // [km^3/s^2] gavitational param. calc'd from Newtons grav. const. & mass of Earth
	m = 2000; // [kg] mass of rocket
	g = 9.81e-3; // [km/s^2]
	Isp = 300; // [s] specific impulse of thruster
	T = 1e4; // thrust [kN]

	// dependent variables
	x = r_t[0];
	y = r_t[1];
	z = r_t[2];
	vx = r_t[3];
	vy = r_t[4];
	vz = r_t[5];
	m = r_t[6];
	r = sqrt(x*x + y*y + z*z);
	v = sqrt(vx*vx + vy*vy + vz*vz);

	// 7 1st order ODEs
	drdt_t[0] = vx;
	drdt_t[1] = vy;
	drdt_t[2] = vz;
	drdt_t[3] = -mu * x / pow(r, 3) + T * vx / (m*v);
	drdt_t[4] = -mu * y / pow(r, 3) + T * vy / (m*v);
	drdt_t[5] = -mu * z / pow(r, 3) + T * vz / (m*v);
	drdt_t[6] = -T / (g*Isp);
}