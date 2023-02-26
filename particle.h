#ifndef PARTICLE_H
#define PARTICLE_H

#include <iostream>

long double rand_ld(long double lower, long double upper);

class Particle											//Base Class of Particle
{
public:
	int id{};
	long double mass{};
	long double x{}, x_prev{}, x_v{}, x_vprev{}, x_accel{}, x_accelprev{};
	long double y{}, y_prev{}, y_v{}, y_vprev{}, y_accel{}, y_accelprev{};
	long double z{}, z_prev{}, z_v{}, z_vprev{}, z_accel{}, z_accelprev{};
	long double density{}, pressure{}, temp_0{ 0 }, temp_f{ 0 };

	Particle(int num, int n)
	{
		if (num != 0) {
			id = num;
			mass = ((long double)5.97e24 / n);
			x = rand_ld(-1, 1) * 5e7,
				y = rand_ld(-1, 1) * 5e7,
				z = rand_ld(-1, 1) * 5e7;
			x_prev = x, y_prev = y, z_prev = z;
			x_vprev = 0, y_vprev = 0, z_vprev = 0;
			int vel = 500;
			if (x > 0) {
				y_vprev = -vel;
			}
			else {
				y_vprev = vel;
			}
			if (y > 0) {
				x_vprev = vel;
			}
			else {
				x_vprev = -vel;
			}
		}
		else {
			id = num;
			mass = (5.97e24 / n);
			x = 0,
				y = 0,
				z = 0;
			x_prev = x, y_prev = y, z_prev = z;
			x_vprev = 0, y_vprev = 0, z_vprev = 0;
			temp_f = 250;
		}
	}
};
#endif