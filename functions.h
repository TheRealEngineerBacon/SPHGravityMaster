#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include "particle.h"
#include <SFML/Graphics.hpp>


//Simulation functions
long double W(long double dist, long double h);

long double grad_W(long double dist, long double h);

void update_pos(Particle** part_array, int n, float delta_t);

std::array<long double, 6> find_centermass(Particle** array, int n, long double h);

void update_grav(Particle** array, int n, short tick, bool friction);

void update_fluid(Particle** array, int n, short tick);

void update_vel(Particle** part_array, int n, float delta_t);

void update_temp(Particle** array, int n, float delta_t, std::array<long double, 6> center_pos);


//Utility functions
void print_Data(Particle** part_array, int n, int focus);

void check_Singularity(Particle** array, int n);

long double rand_ld(long double lower, long double upper);


//Render functions
void update_vertex_pos(Particle** part_array, sf::CircleShape** vertex_array, int n, int pixel_num, 
	float display_radius, float alpha, 
	float beta, float gamma, short tick, 
	int focus, float x_mov, float y_mov,
	std::array<long double, 6> center_pos);

#endif