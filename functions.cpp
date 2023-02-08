#include "functions.h"
#include "particle.h"
#include <SFML/Graphics.hpp>
#include <cstdlib>
#include <cmath>
#include <numeric>
#include <random>
#include <iomanip>
#include <functional>

constexpr long double G = 6.67430e-11;
constexpr long double epsilon = 1e6;
const float sqrt_3 = sqrt(3.f);
const float sqrt_6 = sqrt(6.f);

//Random int gen. Inclusive range (lower <= return <= upper).
long double rand_ld(long double lower, long double upper)
{
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<long double> dist(lower, upper);
	return dist(gen);
}

void update_accel(Particle** array, int n)			//input is pointer to array. Array consists of pointers to objects.
{
	int i{}, j{};
	long double* a_x = new long double[n];
	long double* a_y = new long double[n];
	long double* a_z = new long double[n];
	for (i = 0; i < n; ++i)
	{
		array[i]->x_accelprev = array[i]->x_accel;
		array[i]->y_accelprev = array[i]->y_accel;
		array[i]->z_accelprev = array[i]->z_accel;

		long double i_x{ (*array[i]).x }, i_y{ (*array[i]).y }, i_z{ (*array[i]).z };
		for (j = 0; j < n; ++j)
		{
			long double j_x{ (*array[j]).x }, j_y{ (*array[j]).y }, j_z{ (*array[j]).z };
			if (i != j)
			{
				long double dist = sqrt((j_x - i_x) * (j_x - i_x) + (j_y - i_y) * (j_y - i_y) + (j_z - i_z) * (j_z - i_z));
				long double p{};
				if (dist >= epsilon) {
					p = 1;
				}
				else {
					p = dist / epsilon;
				}
				long double total_accel = G * ((*array[j]).mass / (p * (dist * dist) + (1 - p) * (epsilon * epsilon)));
				a_x[j] = total_accel * ((j_x - i_x) / dist);
				a_y[j] = total_accel * ((j_y - i_y) / dist);
				a_z[j] = total_accel * ((j_z - i_z) / dist);
			}
			else
			{
				a_x[j] = 0;
				a_y[j] = 0;
				a_z[j] = 0;
			}
		}
		(*array[i]).x_accel = std::reduce(a_x, a_x + (n), 0.0l);
		(*array[i]).y_accel = std::reduce(a_y, a_y + (n), 0.0l);
		(*array[i]).z_accel = std::reduce(a_z, a_z + (n), 0.0l);

	}
	delete[] a_x;
	delete[] a_y;
	delete[] a_z;
}

void update_vel(Particle** part_array, int n, float delta_t) {
	int i{};
	for (i = 0; i < n; ++i) {
		(*part_array[i]).x_v = (*part_array[i]).x_vprev + 0.5 * ((*part_array[i]).x_accel + (*part_array[i]).x_accelprev) * delta_t;
		(*part_array[i]).y_v = (*part_array[i]).y_vprev + 0.5 * ((*part_array[i]).y_accel + (*part_array[i]).y_accelprev) * delta_t;
		(*part_array[i]).z_v = (*part_array[i]).z_vprev + 0.5 * ((*part_array[i]).z_accel + (*part_array[i]).z_accelprev) * delta_t;
		(*part_array[i]).x_vprev = (*part_array[i]).x_v;
		(*part_array[i]).y_vprev = (*part_array[i]).y_v;
		(*part_array[i]).z_vprev = (*part_array[i]).z_v;
	}
}

void update_pos(Particle** array, int n, float delta_t) {
	int i{};
	for (i = 0; i < n; ++i) {
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

void update_vertex_pos(Particle** part_array, sf::CircleShape** vertex_array, int n, int pixel_num, float display_radius, float alpha, float beta, float gamma) {
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
	for (i = 0; i < n; ++i) {
		x_raw = static_cast<float>((*part_array[i]).x);
		y_raw = static_cast<float>((*part_array[i]).y);
		z_raw = static_cast<float>((*part_array[i]).z);

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

