#include "functions.h"
#include "particle.h"
#include <SFML/Graphics.hpp>
#include <cstdlib>
#include <cmath>
#include <numeric>
#include <random>
#include <iomanip>
#include <functional>
#include <omp.h>

//Mathematical constants
constexpr long double G = 6.67430e-11;
constexpr long double pi = 3.14159265359;
constexpr long double e = 2.71828182845904523536;
const long double pi_onehalf = pow(pi, 1.5);
const float sqrt_3 = sqrt(3.f);
const float sqrt_6 = sqrt(6.f);
constexpr int thread_n = 8;

//SPH constants
const long double epsilon{ 8e5 };
const long double h = 1e6;
const long double eqstconst = 1e9;
const long double poly_index = 1;

//Random long double generator.
long double rand_ld(long double lower, long double upper)
{
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<long double> dist(lower, upper);
	return dist(gen);
}

long double W(long double dist, long double h) {
	return (1.0 / (h * h * h * pi_onehalf)) * exp(-(dist * dist) / (h * h));
}

long double grad_W(long double dist, long double h) {
	return (-2 / (h * h * h * h * h * pi_onehalf)) * exp(-(dist * dist) / (h * h));
}

void update_accel(Particle** array, int n, short tick)
{
	#pragma omp parallel for num_threads(thread_n)
	for (int i = 0; i < n; ++i)
	{
		long double den{}, pres{}, den_T{}, pres_T{};
		long double a_x{}, a_y{}, a_z{};
		array[i]->x_accelprev = array[i]->x_accel;
		array[i]->y_accelprev = array[i]->y_accel;
		array[i]->z_accelprev = array[i]->z_accel;

		long double i_x{ (*array[i]).x }, i_y{ (*array[i]).y }, i_z{ (*array[i]).z };
		for (int j = 0; j < n; ++j)
		{
			long double j_x{ (*array[j]).x }, j_y{ (*array[j]).y }, j_z{ (*array[j]).z };
			if (i != j)
			{
				long double dist = sqrt((j_x - i_x) * (j_x - i_x) + (j_y - i_y) * (j_y - i_y) + (j_z - i_z) * (j_z - i_z));
				//if (i== 0 && tick % 10000 == 0) {
				//	std::cout << dist << '\n';
				//}
				den = (*array[j]).mass * W(dist, h);
				den_T += den;
				//if (i == 0) {
				//	std::cout << (*array[i]).density << '\n';
				//}
				pres = eqstconst * pow(den, (1 + (1 / poly_index)));
				pres_T += pres;
				//if (i == 0) {
				//	std::cout << "Pressure on P0: " << pressure << '\n';
				//}

				float t{};
				if (dist <= epsilon) {
					t = 1.0f;
				}
				long double total_accel = G * ((*array[j]).mass / ((dist * dist * (1 - t)) + (t * epsilon * epsilon)));
				a_x += total_accel * ((j_x - i_x) / dist);
				a_y += total_accel * ((j_y - i_y) / dist);
				a_z += total_accel * ((j_z - i_z) / dist);
			}
			else
			{
				long double dist = sqrt((j_x - i_x) * (j_x - i_x) + (j_y - i_y) * (j_y - i_y) + (j_z - i_z) * (j_z - i_z));
				den = (*array[j]).mass * (1.0 / (h * h * h * pi_onehalf));
				den_T += den;
				pres = eqstconst * pow(den, (1 + 1 / poly_index));
				pres_T += pres;
			}
		}
		(*array[i]).density = den_T;
		(*array[i]).pressure = pres_T;

		long double nu{ 0.005 };
		(*array[i]).x_accel = a_x - ((*array[i]).x_vprev * nu);
		(*array[i]).y_accel = a_y - ((*array[i]).y_vprev * nu);
		(*array[i]).z_accel = a_z - ((*array[i]).z_vprev * nu);
		//if (i == 0 && tick % 10000 == 0) {
		//	std::cout << "Gravity: " << (*array[i]).x_accel << ", " << (*array[i]).y_accel << ", " << (*array[i]).z_accel << '\n';
		//}
	}

	#pragma omp parallel for num_threads(thread_n)
	for (int i = 0; i < n; ++i) 
	{
		long double fluid_x{}, fluid_y{}, fluid_z{};

		long double i_x{ (*array[i]).x }, i_y{ (*array[i]).y }, i_z{ (*array[i]).z };
		for (int j = 0; j < n; ++j) {
			if (j != i) {
				long double j_x{ (*array[j]).x }, j_y{ (*array[j]).y }, j_z{ (*array[j]).z };
				long double dist = sqrt((j_x - i_x) * (j_x - i_x) + (j_y - i_y) * (j_y - i_y) + (j_z - i_z) * (j_z - i_z));
				
				long double grad_T = grad_W(dist, h);
				long double grad_x = ((j_x - i_x) / dist) * grad_T;
				long double grad_y = ((j_y - i_y) / dist) * grad_T;
				long double grad_z = ((j_z - i_z) / dist) * grad_T;

				long double fluid_force = ((*array[i]).pressure / ((*array[i]).density * (*array[i]).density)) + ((*array[j]).pressure / ((*array[j]).density * (*array[j]).density));
				fluid_x += (*array[j]).mass * fluid_force * grad_x;
				fluid_y += (*array[j]).mass * fluid_force * grad_y;
				fluid_z += (*array[j]).mass * fluid_force * grad_z;

				(*array[i]).x_accel += fluid_x;
				(*array[i]).y_accel += fluid_y;
				(*array[i]).z_accel += fluid_z;
			}
		}
	}
}

void update_vel(Particle** part_array, int n, float delta_t) {
	#pragma omp parallel for num_threads(thread_n)
	for (int i = 0; i < n; ++i) {
		(*part_array[i]).x_v = (*part_array[i]).x_vprev + 0.5 * ((*part_array[i]).x_accel + (*part_array[i]).x_accelprev) * delta_t;
		(*part_array[i]).y_v = (*part_array[i]).y_vprev + 0.5 * ((*part_array[i]).y_accel + (*part_array[i]).y_accelprev) * delta_t;
		(*part_array[i]).z_v = (*part_array[i]).z_vprev + 0.5 * ((*part_array[i]).z_accel + (*part_array[i]).z_accelprev) * delta_t;
		(*part_array[i]).x_vprev = (*part_array[i]).x_v;
		(*part_array[i]).y_vprev = (*part_array[i]).y_v;
		(*part_array[i]).z_vprev = (*part_array[i]).z_v;
	}
}

void update_pos(Particle** array, int n, float delta_t) {
	#pragma omp parallel for num_threads(thread_n)
	for (int i = 0; i < n; ++i) {
		(*array[i]).x = (*array[i]).x_prev + ((*array[i]).x_v * delta_t) + (0.5 * (*array[i]).x_accel * delta_t * delta_t);
		(*array[i]).y = (*array[i]).y_prev + ((*array[i]).y_v * delta_t) + (0.5 * (*array[i]).y_accel * delta_t * delta_t);
		(*array[i]).z = (*array[i]).z_prev + ((*array[i]).z_v * delta_t) + (0.5 * (*array[i]).z_accel * delta_t * delta_t);
		(*array[i]).x_prev = (*array[i]).x;
		(*array[i]).y_prev = (*array[i]).y;
		(*array[i]).z_prev = (*array[i]).z;
	}
}

void print_Energy(Particle** array, int n) {
	int i{}, j{};
	int array_size = n;
	long double E_pot{};
	long double E_kin{};
	long double E_total{};
	for (i = 0; i < array_size; ++i)
	{
		long double vel = sqrt(pow(abs((*array[i]).x_v), 2) + pow(abs((*array[i]).y_v), 2) + pow(abs((*array[i]).z_v), 2));
		E_kin += 0.5 * (*array[i]).mass * pow(vel, 2);

		long double i_x{ (*array[i]).x }, i_y{ (*array[i]).y }, i_z{ (*array[i]).z };
		for (j = i + 1; j < array_size; ++j)
		{
			long double j_x{ (*array[j]).x }, j_y{ (*array[j]).y }, j_z{ (*array[j]).z };
			{
				long double dist = sqrt(pow(abs(j_x - i_x), 2) + pow(abs(j_y - i_y), 2) + pow(abs(j_z - i_z), 2));
				E_pot += -G * (((*array[i]).mass * (*array[j]).mass) / dist);
			}
		}
	}
	E_total = E_pot + E_kin;
	std::cout << std::setprecision(16) << E_total << '\n';
}

void update_vertex_pos(Particle** array, sf::CircleShape** vertex_array, int n, int pixel_num, float display_radius, float alpha, float beta, float gamma, short tick, int focus) {
	int i{};
	float x_raw{}, y_raw{}, z_raw{};
	float x_pos{}, y_pos{}, x_adj{}, y_adj{};
	const float sin_a = std::sin(alpha);
	const float sin_b = std::sin(beta);
	const float sin_g = std::sin(gamma);
	const float cos_a = std::cos(alpha);
	const float cos_b = std::cos(beta);
	const float cos_g = std::cos(gamma);
	const float radius = (*vertex_array[0]).getRadius();
	const float scale = display_radius / pixel_num;
	const float middle = static_cast<float>(pixel_num) / 2;

	float cenx, ceny, cenz;
	cenx = static_cast<float>((*array[focus]).x);
	ceny = static_cast<float>((*array[focus]).y);
	cenz = static_cast<float>((*array[focus]).z);

	for (i = 0; i < n; ++i) {
		x_raw = static_cast<float>((*array[i]).x) - cenx;
		y_raw = static_cast<float>((*array[i]).y) - ceny;
		z_raw = static_cast<float>((*array[i]).z) - cenz;

		//Rotation then orthographic projection.
		x_adj = x_raw * (cos_b * cos_g)
			+ y_raw * ((sin_a * sin_b * cos_g) - (cos_a * sin_g))
			+ z_raw * ((cos_a * sin_b * cos_g) + (sin_a * sin_g));

		y_adj = x_raw * (cos_b * sin_g)
			+ y_raw * ((sin_a * sin_b * sin_g) + (cos_a * cos_g))
			+ z_raw * ((cos_a * sin_b * sin_g) - (sin_a * cos_g));

		x_pos = (x_adj / scale) + middle - radius;
		y_pos = (y_adj / scale) + middle - radius;

		(*vertex_array[i]).setPosition(x_pos, y_pos);
	}
}

