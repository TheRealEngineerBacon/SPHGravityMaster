#include "functions.h"
#include "particle.h"
#include <SFML/Graphics.hpp>
#include <array>
#include <random>
#include <iomanip>
#include <omp.h>
#include "SmoothGrid.h"

#define MULTI

//Mathematical constants
constexpr long double G = -6.67430e-11;
constexpr long double pi = 3.14159265359;
constexpr long double e = 2.718281828459045;
const long double one_over_pi = 1 / pi;
constexpr int thread_n = 20;

//SPH constants
constexpr long double epsilon{ 1.5e6 };
constexpr long double epsilon_squared{ epsilon * epsilon };
constexpr long double h = 1e6;
constexpr long double eqstconst = 0.1;


//Simulation functions
void set_initial_position(Particle** array, int n) {
	int max = (int)ceil(cbrt(n));
	int mid = (int)floor((max) / 2);
	long double length = 2.5e6;
	for (int p = 0; p < n; ++p) {
		int z = p % max;
		int y = ((p - z) / max) % max;
		int x = (((p - z) / max) - y) / max;
		(*array[p]).x = (x - mid) * length;
		(*array[p]).y = (y - mid) * length;
		(*array[p]).z = (z - mid) * length;
		(*array[p]).x_prev = (*array[p]).x;
		(*array[p]).y_prev = (*array[p]).y;
		(*array[p]).z_prev = (*array[p]).z;
	}
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
	else {
		std::cout << "Error: Function W: Negative q value encountered." << '\n';
		return 0.0L;
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
		const long double sigma = 0.75 * one_over_pi * h1 * h1 * h1;
		const long double inner = 2.0 - q;
		return sigma * inner * inner;
	}
	else if (q >= 0) {
		const long double sigma = one_over_pi * h1 * h1 * h1;
		return sigma * (-3 * q * (1 - 0.75 * q));
	}
	else {
		std::cout << "Error: Function grad_W: Negative q value encountered." << '\n';
		return 0.0L;
	}
	//return (-2 / (h * h * h * h * h * pi_threehalf)) * exp(-(dist * dist) / (h * h));
}

void update_pos(Particle** array, int n, double delta_t) {
#ifdef MULTI
	#pragma omp parallel for num_threads(thread_n)
#endif
	for (int i = 0; i < n; ++i) {
		(*array[i]).x = (*array[i]).x_prev + ((*array[i]).x_v * delta_t) + (0.5 * (*array[i]).x_accel * delta_t * delta_t);
		(*array[i]).y = (*array[i]).y_prev + ((*array[i]).y_v * delta_t) + (0.5 * (*array[i]).y_accel * delta_t * delta_t);
		(*array[i]).z = (*array[i]).z_prev + ((*array[i]).z_v * delta_t) + (0.5 * (*array[i]).z_accel * delta_t * delta_t);
		(*array[i]).x_prev = (*array[i]).x;
		(*array[i]).y_prev = (*array[i]).y;
		(*array[i]).z_prev = (*array[i]).z;
	}
}

std::array<long double, 6> find_centermass(Particle** array, int n, long double h) {
	std::array<long double, 6> temp;
	long double temp_x{}, temp_y{}, temp_z{}, mass_total{};
	long double x_absmax{}, y_absmax{}, z_absmax{};
	for (int i = 0; i < n; ++i) {
		mass_total += (*array[i]).mass;
		temp_x += (*array[i]).mass * (*array[i]).x;
		temp_y += (*array[i]).mass * (*array[i]).y;
		temp_z += (*array[i]).mass * (*array[i]).z;
	}
	temp[0] = temp_x / mass_total;
	temp[1] = temp_y / mass_total;
	temp[2] = temp_z / mass_total;
	
	for (int i = 0; i < n; ++i) {
		long double x_abs{}, y_abs{}, z_abs{};
		x_abs = std::abs((*array[i]).x - temp[0]);
		y_abs = std::abs((*array[i]).y - temp[1]);
		z_abs = std::abs((*array[i]).z - temp[2]);
		if (x_abs > x_absmax)
			x_absmax = x_abs;
		if (y_abs > y_absmax)
			y_absmax = y_abs;
		if (z_abs > z_absmax)
			z_absmax = z_abs;
	}
	temp[3] = x_absmax + 1e3;
	temp[4] = y_absmax + 1e3;
	temp[5] = z_absmax + 1e3;
	return temp;
}

void update_surr(Particle** array, int n, int grid_size, SmoothGrid* grid) {
	#pragma omp parallel for num_threads(thread_n)
	for (int p = 0; p < n; ++p) {
		(*array[p]).surr.clear();
		const int p_x{ (*array[p]).ind_x }, p_y{ (*array[p]).ind_y }, p_z{ (*array[p]).ind_z };
		
		for (int i = p_x - 1; i < (p_x + 2); ++i) {
			if (i < 0 || i > grid_size - 1) {
				continue;
			}

			for (int j = p_y - 1; j < (p_y + 2); ++j) {
				if (j < 0 || j > grid_size - 1) {
					continue;
				}

				for (int k = p_z - 1; k < (p_z + 2); ++k) {
					if (k < 0 || k > grid_size - 1) {
						continue;
					}
					if ((*grid).grid[i][j][k].empty()) {
						continue;
					}
					for (int m = 0; m < (*grid).grid[i][j][k].size(); ++m)
						(*array[p]).surr.emplace_back((*grid).grid[i][j][k][m]);
				}
			}
		}
	}
}

void update_grav(Particle** array, int n, short tick, bool friction) {
#ifdef MULTI
	#pragma omp parallel for num_threads(thread_n)
#endif
	for (int i = 0; i < n; ++i)
	{
		long double a_x{}, a_y{}, a_z{};
		array[i]->x_accelprev = array[i]->x_accel;
		array[i]->y_accelprev = array[i]->y_accel;
		array[i]->z_accelprev = array[i]->z_accel;

		const long double i_x{ (*array[i]).x }, i_y{ (*array[i]).y }, i_z{ (*array[i]).z };

		for (int j = 0; j < n; ++j)
		{
			const long double j_x{ (*array[j]).x }, j_y{ (*array[j]).y }, j_z{ (*array[j]).z };
			if (i != j)
			{
				const long double diff = (j_x - i_x) * (j_x - i_x) + (j_y - i_y) * (j_y - i_y) + (j_z - i_z) * (j_z - i_z);
				const long double dist = sqrtl(diff);

				int t{};
				if (dist <= epsilon) {
					t = 1;
				}
				else {
					t = 0;
				}
				
				const long double total_accel = G * (*array[j]).mass * (1 / (((diff * (1 - t)) + (epsilon_squared * t)) * dist));
				a_x += total_accel * (i_x - j_x);
				a_y += total_accel * (i_y - j_y);
				a_z += total_accel * (i_z - j_z);
			}
		}

		long double nu_x{}, nu_y{}, nu_z{};
		if (friction == true) {
			nu_x = 0.001; 
			nu_y = 0.001; 
			nu_z = 0.001;
		}
		else {
			nu_x = 0.0; 
			nu_y = 0.0; 
			nu_z = 0.0;
		}
		const long double b = 1e-2;
		(*array[i]).x_accel = a_x - ((*array[i]).x_vprev * nu_x) + (-b * ((*array[i]).temp_f - 130) * a_x);
		(*array[i]).y_accel = a_y - ((*array[i]).y_vprev * nu_y) + (-b * ((*array[i]).temp_f - 130) * a_y);
		(*array[i]).z_accel = a_z - ((*array[i]).z_vprev * nu_z) + (-b * ((*array[i]).temp_f - 130) * a_z);
	}
}

void update_density(Particle** array, int n) {
#ifdef MULTI
	#pragma omp parallel for num_threads(thread_n)
#endif
	for (int i = 0; i < n; ++i)
	{
		const long double i_x{ (*array[i]).x }, i_y{ (*array[i]).y }, i_z{ (*array[i]).z };
		long double den_T{};

		for (int j : (*array[i]).surr)
		{
			if (j != i) {
				const long double j_x{ (*array[j]).x }, j_y{ (*array[j]).y }, j_z{ (*array[j]).z };
				const long double dist = sqrt((j_x - i_x) * (j_x - i_x) + (j_y - i_y) * (j_y - i_y) + (j_z - i_z) * (j_z - i_z));

				den_T += (*array[j]).mass * W(dist, h);

			}
			else {
				den_T += (*array[j]).mass * (1.0 / (h * h * h * pi));
			}
		}
		(*array[i]).density = den_T;
	}
}

void update_fluid(Particle** array, int n, short tick) {
	//Computation of fluid pressures from density values.
#ifdef MULTI
	#pragma omp parallel for num_threads(thread_n)
#endif
	for (int i = 0; i < n; ++i) {
		(*array[i]).pressure = eqstconst * (*array[i]).density * (*array[i]).density;
	}

	//Calculate acceleration due to fluid pressure.
#ifdef MULTI
	#pragma omp parallel for num_threads(thread_n)
#endif
	for (int i = 0; i < n; ++i)
	{
		long double fluid_x{}, fluid_y{}, fluid_z{};
		const long double i_x{ (*array[i]).x }, i_y{ (*array[i]).y }, i_z{ (*array[i]).z };

		for (int j : (*array[i]).surr) {
			if (j != i) {
				const long double j_x{ (*array[j]).x }, j_y{ (*array[j]).y }, j_z{ (*array[j]).z };
				const long double dist = sqrt((j_x - i_x) * (j_x - i_x) + (j_y - i_y) * (j_y - i_y) + (j_z - i_z) * (j_z - i_z));

				const long double grad_T = grad_W(dist, h);
				const long double grad_x = ((i_x - j_x) / dist) * grad_T;
				const long double grad_y = ((i_y - j_y) / dist) * grad_T;
				const long double grad_z = ((i_z - j_z) / dist) * grad_T;

				const long double fluid_force = ((*array[i]).pressure / ((*array[i]).density * (*array[i]).density)) + ((*array[j]).pressure / ((*array[j]).density * (*array[j]).density));
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

void update_vel(Particle** part_array, int n, double delta_t) {
#ifdef MULTI
	#pragma omp parallel for num_threads(thread_n)
#endif
	for (int i = 0; i < n; ++i) {
		(*part_array[i]).x_v = (*part_array[i]).x_vprev + 0.5 * ((*part_array[i]).x_accel + (*part_array[i]).x_accelprev) * delta_t;
		(*part_array[i]).y_v = (*part_array[i]).y_vprev + 0.5 * ((*part_array[i]).y_accel + (*part_array[i]).y_accelprev) * delta_t;
		(*part_array[i]).z_v = (*part_array[i]).z_vprev + 0.5 * ((*part_array[i]).z_accel + (*part_array[i]).z_accelprev) * delta_t;
		(*part_array[i]).x_vprev = (*part_array[i]).x_v;
		(*part_array[i]).y_vprev = (*part_array[i]).y_v;
		(*part_array[i]).z_vprev = (*part_array[i]).z_v;
	}
}

void update_temp(Particle** array, int n, double delta_t, std::array<long double, 6> center_pos) {
#ifdef MULTI
	#pragma omp parallel for num_threads(thread_n)
#endif
	for (int i = 0; i < n; ++i) {
		(*array[i]).temp_0 = (*array[i]).temp_f;
		const long double i_x{ (*array[i]).x }, i_y{ (*array[i]).y }, i_z{ (*array[i]).z };
		
		const long double dist_from_center = sqrt((center_pos[0] - i_x) * (center_pos[0] - i_x)
			+ (center_pos[1] - i_y) * (center_pos[1] - i_y)
			+ (center_pos[2] - i_z) * (center_pos[2] - i_z));
		long double delta_temp{};
		
		for (int j : (*array[i]).surr)
		{
			const long double j_x{ (*array[j]).x }, j_y{ (*array[j]).y }, j_z{ (*array[j]).z };
			if (j != i)
			{
				const long double dist = sqrt((j_x - i_x) * (j_x - i_x) + (j_y - i_y) * (j_y - i_y) + (j_z - i_z) * (j_z - i_z));
				const long double temp_diff{ (*array[j]).temp_f - (*array[i]).temp_f };
				const long double gradient{ grad_W(dist, h) };
				const long double mass_density{ (*array[j]).mass / (*array[j]).density };
				delta_temp += mass_density * temp_diff * gradient;
			}
		}
		const long double conduc{ 1000 };
		const long double c_p{ 1 };
		if (dist_from_center > (center_pos[3] - (0.5 * h))) {
			(*array[i]).temp_f  = 30;
		} else if (dist_from_center > (2 * h)) {
			(*array[i]).temp_f += ((2 * conduc)/((*array[i]).density * c_p)) * delta_temp * delta_t;
		} else if (dist_from_center <= (2 * h)) {
			(*array[i]).temp_f = 240;
		}
	}
}


//Utility functions
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

long double rand_ld(long double lower, long double upper)
{
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<long double> dist(lower, upper);
	return dist(gen);
}


//Render functions
void update_vertex_pos(Particle** array, sf::CircleShape** vertex_array, int n, int pixel_num, 
	float display_radius, float alpha, 
	float beta, float gamma, short tick, 
	int focus, float x_mov, float y_mov,
	std::array<long double, 6> center_pos) {
	
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

	long double max_den{};
	for (int i = 0; i < n; ++i) {
		if ((*array[i]).density > max_den) {
			max_den = (*array[i]).density;
		}

		x_raw = static_cast<float>((*array[i]).x - center_pos[0]);
		y_raw = static_cast<float>((*array[i]).y - center_pos[1]);
		z_raw = static_cast<float>((*array[i]).z - center_pos[2]);

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
			(*vertex_array[i]).setRadius(static_cast<float>(((*array[i]).density / (max_den + 0.01)) * 2));
		}
		else {
			(*vertex_array[i]).setRadius(static_cast<float>(((*array[i]).density / (max_den + 0.01)) * 4));
		}
		
		//long double value = abs(((*array[i]).density - min_den) / (max_den - min_den));
		sf::Color color{ static_cast<uint8_t>((*array[i]).temp_0), 0, static_cast<uint8_t>(255 - (*array[i]).temp_0), 255 };
		(*vertex_array[i]).setFillColor(color);
	}
}





