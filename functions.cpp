#include "functions.h"
#include "particle.h"
#include <SFML/Graphics.hpp>
#include <random>
#include <iomanip>
#include <omp.h>

//Mathematical constants
constexpr long double G = 6.67430e-11;
constexpr long double pi = 3.14159265359;
constexpr long double e = 2.718281828459045;
const long double one_over_pi = 1 / pi;
constexpr int thread_n = 16;

//SPH constants
constexpr long double epsilon{ 1.4e6 };
constexpr long double h = 1e6;
constexpr long double eqstconst = 0.1;

//Random long double generator.
long double rand_ld(long double lower, long double upper)
{
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<long double> dist(lower, upper);
	return dist(gen);
}

long double W(long double dist, long double h) {
	const long double h1 = 1.0L / h;
	const long double q = dist * h1;
	if (q > 2) {
		return 0.0L;
	}
	else if (q >= 1) {
		const long double sigma = 0.25 * one_over_pi * h1 * h1 * h1;
		const long double inner = 2.0 - q;
		return sigma * inner * inner * inner;
	}
	else if (q >= 0) {
		const long double sigma = one_over_pi * h1 * h1 * h1;
		return sigma * (1 - (1.5 * q * q) + (0.75 * q * q * q));
	}
	//return (1.0 / (h * h * h * pi_threehalf)) * exp(-(dist * dist) / (h * h));
}

long double grad_W(long double dist, long double h) {
	const long double h1 = 1.0L / h;
	const long double q = dist * h1;
	if (q > 2) {
		return 0.0L;
	}
	else if (q >= 1) {
		long double sigma = 0.75 * one_over_pi * h1 * h1 * h1;
		long double inner = 2.0 - q;
		return sigma * inner * inner;
	}
	else if (q >= 0) {
		long double sigma = one_over_pi * h1 * h1 * h1;
		return sigma * (-3 * q * (1 - 0.75 * q));
	}
	//return (-2 / (h * h * h * h * h * pi_threehalf)) * exp(-(dist * dist) / (h * h));
}

void update_temp(Particle** array, int n, float delta_t) {
	#pragma omp parallel for num_threads(thread_n)
	for (int i = 0; i < n; ++i) {
		(*array[i]).temp_0 = (*array[i]).temp_f;
		long double delta_temp{};
		
		long double i_x{ (*array[i]).x }, i_y{ (*array[i]).y }, i_z{ (*array[i]).z };
		for (int j = 0; j < n; ++j)
		{
			long double j_x{ (*array[j]).x }, j_y{ (*array[j]).y }, j_z{ (*array[j]).z };
			if (j != i)
			{
				long double dist = sqrt((j_x - i_x) * (j_x - i_x) + (j_y - i_y) * (j_y - i_y) + (j_z - i_z) * (j_z - i_z));
				long double temp_diff{ (*array[j]).temp_f - (*array[i]).temp_f };
				long double gradient{ grad_W(dist, h) };
				long double mass_density{ (*array[j]).mass / (*array[j]).density };
				delta_temp += mass_density * temp_diff * gradient;
			}
		}
		long double conduc{ 1000 };
		long double c_p{ 10 };
		if (i != 0) {
			(*array[i]).temp_f += ((2 * conduc)/((*array[i]).density * c_p)) * delta_temp * delta_t;
		}
		else {
			(*array[i]).temp_f = 240;
		}
	}
}

void update_grav(Particle** array, int n, short tick)
{
	//Computation of gravitation forces and density summation.
	#pragma omp parallel for num_threads(thread_n)
	for (int i = 0; i < n; ++i)
	{
		long double a_x{}, a_y{}, a_z{}, den_T{};
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

				den_T += (*array[j]).mass * W(dist, h);

				float t{};
				if (dist <= epsilon) {
					t = 1.0f;
				}
				else {
					t = 0.0f;
				}
				long double total_accel = -G * ((*array[j]).mass / ((dist * dist * (1 - t)) + (t * epsilon * epsilon)));
				a_x += total_accel * ((i_x - j_x) / dist);
				a_y += total_accel * ((i_y - j_y) / dist);
				a_z += total_accel * ((i_z - j_z) / dist);

			}
			else
			{
				den_T += (*array[j]).mass * (1.0 / (h * h * h * pi));
			}
		}
		(*array[i]).density = den_T;

		long double nu_x{ 0.0002 };
		long double nu_y{ 0.0002 };
		long double nu_z{ 0.0002 };
		(*array[i]).x_accel = a_x - ((*array[i]).x_vprev * nu_x);
		(*array[i]).y_accel = a_y - ((*array[i]).y_vprev * nu_y);
		(*array[i]).z_accel = a_z - ((*array[i]).z_vprev * nu_z);
	}
}

void update_fluid(Particle** array, int n, short tick) {
	//Computation of fluid pressures from density values.
	#pragma omp parallel for num_threads(thread_n)
	for (int i = 0; i < n; ++i) {
		(*array[i]).pressure = eqstconst * (*array[i]).density * (*array[i]).density;
	}

	//Calculate acceleration due to fluid pressure.
	#pragma omp parallel for num_threads(thread_n)
	for (int i = 0; i < n; ++i)
	{
		long double fluid_x{}, fluid_y{}, fluid_z{};
		long double i_x{ (*array[i]).x }, i_y{ (*array[i]).y }, i_z{ (*array[i]).z };

		for (int j = 0; j < n; ++j) {
			if (j != i) {
				long double j_press = eqstconst * (*array[j]).density * (*array[j]).density;
				long double j_x{ (*array[j]).x }, j_y{ (*array[j]).y }, j_z{ (*array[j]).z };
				long double dist = sqrt((j_x - i_x) * (j_x - i_x) + (j_y - i_y) * (j_y - i_y) + (j_z - i_z) * (j_z - i_z));

				long double grad_T = grad_W(dist, h);
				long double grad_x = ((i_x - j_x) / dist) * grad_T;
				long double grad_y = ((i_y - j_y) / dist) * grad_T;
				long double grad_z = ((i_z - j_z) / dist) * grad_T;

				long double fluid_force = ((*array[i]).pressure / ((*array[i]).density * (*array[i]).density)) + ((*array[j]).pressure / ((*array[j]).density * (*array[j]).density));
				fluid_x += (*array[j]).mass * fluid_force * grad_x;
				fluid_y += (*array[j]).mass * fluid_force * grad_y;
				fluid_z += (*array[j]).mass * fluid_force * grad_z;

			}
		}
		(*array[i]).x_accel += fluid_x;
		(*array[i]).y_accel += fluid_y;
		(*array[i]).z_accel += fluid_z;
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

void print_Data(Particle** array, int n, int focus) {
	long double E_pot{};
	long double E_kin{};
	long double E_total{};
	//for (int i = 0; i < n; ++i)
	//{
	//	long double vel = sqrt(pow(abs((*array[i]).x_v), 2) + pow(abs((*array[i]).y_v), 2) + pow(abs((*array[i]).z_v), 2));
	//	E_kin += 0.5 * (*array[i]).mass * pow(vel, 2);

	//	long double i_x{ (*array[i]).x }, i_y{ (*array[i]).y }, i_z{ (*array[i]).z };
	//	for (int j = i + 1; j < n; ++j)
	//	{
	//		long double j_x{ (*array[j]).x }, j_y{ (*array[j]).y }, j_z{ (*array[j]).z };
	//		{
	//			long double dist = sqrt(pow(abs(j_x - i_x), 2) + pow(abs(j_y - i_y), 2) + pow(abs(j_z - i_z), 2));
	//			E_pot += -G * (((*array[i]).mass * (*array[j]).mass) / dist);
	//		}
	//	}
	//}
	//E_total = E_pot + E_kin;
	std::cout << std::setprecision(12) << (*array[focus]).density << '\n';
}

void update_vertex_pos(Particle** array, sf::CircleShape** vertex_array, int n, int pixel_num, float display_radius, float alpha, float beta, float gamma, short tick, int focus, float x_mov, float y_mov) {
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

	long double mass_x{}, mass_y{}, mass_z{}, mass_total{};
	for (int i = 0; i < n; ++i) {
		mass_total += (*array[i]).mass;
		mass_x = (*array[i]).mass * (*array[i]).x;
		mass_y = (*array[i]).mass * (*array[i]).y;
		mass_z = (*array[i]).mass * (*array[i]).z;
	}

	//Shifts coordinates to focus on specific particle.
	//float cenx, ceny, cenz;
	//cenx = static_cast<float>((*array[focus]).x);
	//ceny = static_cast<float>((*array[focus]).y);
	//cenz = static_cast<float>((*array[focus]).z);

	//Shifts coordinates to focus on center of mass.
	float cenx, ceny, cenz;
	cenx = static_cast<float>(mass_x / mass_total);
	ceny = static_cast<float>(mass_y / mass_total);
	cenz = static_cast<float>(mass_z / mass_total);

	long double max_den{};
	for (int i = 0; i < n; ++i) {
		if ((*array[i]).density > max_den) {
			max_den = (*array[i]).density;
		}

		x_raw = static_cast<float>((*array[i]).x) - cenx;
		y_raw = static_cast<float>((*array[i]).y) - ceny;
		z_raw = static_cast<float>((*array[i]).z) - cenz;

		//Rotation then orthographic projection.
		x_adj = x_raw * (cos_b * cos_g)
			+ y_raw * ((sin_a * sin_b * cos_g) - (cos_a * sin_g))
			+ z_raw * ((cos_a * sin_b * cos_g) + (sin_a * sin_g))
			+ x_mov;


		y_adj = x_raw * (cos_b * sin_g)
			+ y_raw * ((sin_a * sin_b * sin_g) + (cos_a * cos_g))
			+ z_raw * ((cos_a * sin_b * sin_g) - (sin_a * cos_g))
			+ y_mov;

		x_pos = (x_adj / scale) + middle - radius;
		y_pos = (y_adj / scale) + middle - radius;

		(*vertex_array[i]).setPosition(x_pos, y_pos);
	}

	for (int i = 0; i < n; ++i) {
		if (i != focus) {
			(*vertex_array[i]).setRadius(static_cast<float>(((*array[i]).density / max_den) * 2));
		}
		else {
			(*vertex_array[i]).setRadius(static_cast<float>(((*array[i]).density / max_den) * 4));
		}
		
		//long double value = abs(((*array[i]).density - min_den) / (max_den - min_den));
		sf::Color color{ static_cast<uint8_t>((*array[i]).temp_0), 0, 0, 255 };
		(*vertex_array[i]).setFillColor(color);
	}
}

void check_Singularity(Particle** array, int n) {
	/* 
	This function checks for singularities, events in which the gravitiational forces overcomes the fluid force.
	The result is two particles beings sucked together, thus occupying the same location in space. This is not desirable.
	*/
	bool singularity = false;
	for (int i = 0; i < n; ++i) {
		long double i_x{ (*array[i]).x }, i_y{ (*array[i]).y }, i_z{ (*array[i]).z };

		for (int j = i+1; j < n; ++j) {
			long double j_x{ (*array[j]).x }, j_y{ (*array[j]).y }, j_z{ (*array[j]).z };
			long double dist = sqrt((j_x - i_x) * (j_x - i_x) + (j_y - i_y) * (j_y - i_y) + (j_z - i_z) * (j_z - i_z));
			if (dist < 1e5) {
				singularity = true;
			}
		}
	}
	if (singularity == true) {
		std::cout << "Warning: Singularity event detected." << '\n';
	}
}