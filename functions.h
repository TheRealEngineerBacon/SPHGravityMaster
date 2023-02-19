#ifndef FUNCTIONS_H
#define FUNCTIONS_H
#include "particle.h"
#include <SFML/Graphics.hpp>

long double rand_ld(long double lower, long double upper);

long double W(long double dist, long double h);

void update_accel(Particle** array, int n, short tick);

void update_vel(Particle** part_array, int n, float delta_t);

void update_pos(Particle** part_array, int n, float delta_t);

void print_Data(Particle** part_array, int n, int focus);

void update_vertex_pos(Particle** part_array, sf::CircleShape** vertex_array, int n, int pixel_num, float display_radius, float alpha, float beta, float gamma, short tick, int focus);

void check_Singularity(Particle** array, int n);
#endif