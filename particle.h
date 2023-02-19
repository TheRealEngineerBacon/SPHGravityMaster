#ifndef PARTICLE_H
#define PARTICLE_H

#include "particle.h"									//#include and forwards here.
#include <iostream>
#include <sstream>
#include <fstream>

long double rand_ld(long double lower, long double upper);

class Particle											//Base Class of Particle
{
public:
	int id{};
	long double mass{};
	long double x{}, x_prev{}, x_v{}, x_vprev{}, x_accel{}, x_accelprev{};
	long double y{}, y_prev{}, y_v{}, y_vprev{}, y_accel{}, y_accelprev{};
	long double z{}, z_prev{}, z_v{}, z_vprev{}, z_accel{}, z_accelprev{};
	long double density{}, pressure{};

	Particle(int num, int n)
	{
		if (num != 0) {
			id = num;
			mass = (5.97e24 / n) + (rand_ld(-1, 1) * 2e22);
			x = rand_ld(-1, 1) * 3e7,
				y = rand_ld(-1, 1) * 3e7,
				z = rand_ld(-1, 1) * 3e7;
			//3e7
			x_prev = x, y_prev = y, z_prev = z;
			x_vprev = 0, y_vprev = 0, z_vprev = 0;
			//rand_ld(-1, 1) * 6.371e6
		}
		else {
			id = num;
			mass = (5.97e24 / n) + (rand_ld(-1, 1) * 2e22);
			x = 0,
				y = 0,
				z = 0;
			x_prev = x, y_prev = y, z_prev = z;
			x_vprev = 0, y_vprev = 0, z_vprev = 0;
		}
		

		//if (num == 0) {
		//	id = num;
		//	mass = 5.97e24 / n;
		//	x = 3e6,
		//		y = 0,
		//		z = 0;
		//	x_prev = x, y_prev = y, z_prev = z;
		//	x_vprev = 0, y_vprev = 0, z_vprev = 0;
		//}
		//else {
		//	id = num;
		//	mass = 5.97e24 / n;
		//	x = -3e6,
		//		y = 0,
		//		z = 0;
		//	x_prev = x, y_prev = y, z_prev = z;
		//	x_vprev = 0, y_vprev = 0, z_vprev = 0;
		//}
	}
};
#endif