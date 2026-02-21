#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <nummethods.h>
#include <math.h>

void orbitalburn(double drdt_t[], double r_t[], double t, uint Nr) {
	(void) t; // time independent
	if (Nr != 7) {
		fprintf(stderr, "ERROR: Given %d dependent var values, should be 2.\n", Nr);
		exit(1);
	}
	// constants
	double mu = 3.986e8; // [km^3/s^2] gavitational param. calc'd from Newtons grav. const. & mass of Earth
	double g = 9.81e-3; // [km/s^2]
	double Isp = 300; // [s] specific impulse of thruster
	double T = 1e4; // thrust [kN]

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
	double mu = 3.986e8; // [km^3/s^2] gavitational param. calc'd from Newtons grav. const. & mass of Earth
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