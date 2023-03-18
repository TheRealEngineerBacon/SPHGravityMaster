#ifndef CONSTANTS_H
#define CONSTANTS_H
#include <array>

//Particle Constants


//Simulation Constants
constexpr int n{ 512 };
constexpr int compute_ratio = 100;
std::array<long double, 6> center_pos{};
double delta_t = 0.2;
short tick{};
bool friction{ false };
constexpr long double h = 1e6;
const int grid_size = 17;


//Render Constants
constexpr int pixel_num = 512;
float alpha{}, beta{}, gamma{}, scale{};
int focus{};
float center_mass{}, x_mov, y_mov;
float display_radius{ 1e8 };
bool pause_state = false;

#endif