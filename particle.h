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

	//long double x_hist[101]{};						//Stores history for file output/trailing.
	//long double y_hist[101]{};
	//long double z_hist[101]{};

	Particle(int num)
	{
		id = num;
		mass = 5.97e24;
		x = rand_ld(-1, 1) * 3.5e8,
			y = rand_ld(-1, 1) * 3.5e8,
			z = rand_ld(-1, 1) * 3.5e8;
		x_prev = x, y_prev = y, z_prev = z;
		x_vprev = 0, y_vprev = 0, z_vprev = 0;
	}

	//x_hist[0] = x;
	//y_hist[0] = y;
	//z_hist[0] = z;
};

#endif