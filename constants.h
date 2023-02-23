#ifndef CONSTANTS_H
#define CONSTANTS_H

//Particle Constants


//Simulation Constants
constexpr int n{ 160 };
constexpr float delta_t = 0.5;
constexpr int compute_ratio = 200;
short tick{};


//Render Constants
constexpr int pixel_num = 512;
float alpha{}, beta{}, gamma{}, scale{};
int focus{};
float center_mass{}, x_mov, y_mov;
float display_radius{ 1e8 };

#endif