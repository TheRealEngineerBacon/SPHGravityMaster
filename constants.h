#ifndef CONSTANTS_H
#define CONSTANTS_H
#include <array>

//Particle Constants


//Simulation Constants
constexpr int n{ 192 };
constexpr int compute_ratio = 200;
std::array<long double, 3> center_pos{};
float delta_t = 0.5;
short tick{};
bool friction{ false };


//Render Constants
constexpr int pixel_num = 768;
float alpha{}, beta{}, gamma{}, scale{};
int focus{};
float center_mass{}, x_mov, y_mov;
float display_radius{ 1e8 };
bool pause_state = false;

#endif