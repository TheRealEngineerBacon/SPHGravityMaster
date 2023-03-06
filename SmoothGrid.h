#ifndef SMOOTHGRID_H
#define SMOOTHGRID_H

#include "particle.h"
#include <vector>
#include <array>

class SmoothGrid {
public:
	std::vector<std::vector<std::vector<std::vector<int>>>> grid;

	SmoothGrid(int size, long double h, int n) {
		grid.resize(size, std::vector<std::vector<std::vector<int>>>
				   (size, std::vector<std::vector<int>>
				   (size, std::vector<int>())));
	}

	void update_grid(Particle** array, const int n, const int grid_size, const long double h, const std::array<long double, 6> center_pos) {
		long double h_adj = calculate_h(grid_size, h, center_pos);
		for (int i = 0; i < n; ++i) {
			if ((*array[i]).on_grid == true) {
				const int offset = (int)floor(grid_size / 2.f);
				const long double x{ (*array[i]).x }, y{ (*array[i]).y }, z{ (*array[i]).z };
				const long double h_half = h_adj / 2;
				const int a = (int)floor((x + h_half - center_pos[0]) / h_adj) + offset;
				const int b = (int)floor((y + h_half - center_pos[1]) / h_adj) + offset;
				const int c = (int)floor((z + h_half - center_pos[2]) / h_adj) + offset;
				if ((*array[i]).ind_x != a || (*array[i]).ind_y != b || (*array[i]).ind_z != c) {
					int ind_x = (*array[i]).ind_x;
					int ind_y = (*array[i]).ind_y;
					int ind_z = (*array[i]).ind_z;
					int id = (*array[i]).id;
					int erasure = find_element(id, ind_x, ind_y, ind_z);
					grid[ind_x][ind_y][ind_z].erase(grid[ind_x][ind_y][ind_z].begin() + erasure);
					grid[a][b][c].push_back((*array[i]).id);
					(*array[i]).ind_x = a;
					(*array[i]).ind_y = b;
					(*array[i]).ind_z = c;
				}
			}
			else {
				place(array, i, grid_size, h, center_pos, h_adj);
				(*array[i]).on_grid = true;
			}
		}
	}

	void place(Particle** array, const int i, const int grid_size, const long double h, const std::array<long double, 6> center_pos, const long double h_adj) {
		const int offset = (int)floor(grid_size / 2.f);
		const long double x{ (*array[i]).x }, y{ (*array[i]).y }, z{ (*array[i]).z };
		const long double h_half = h_adj / 2;
		const int a = (int)floor((x + h_half - center_pos[0]) / h_adj) + offset;
		const int b = (int)floor((y + h_half - center_pos[1]) / h_adj) + offset;
		const int c = (int)floor((z + h_half - center_pos[2]) / h_adj) + offset;
		grid[a][b][c].push_back((*array[i]).id);
		(*array[i]).ind_x = a;
		(*array[i]).ind_y = b;
		(*array[i]).ind_z = c;
	}

	int find_element(const int id, const int ind_x, const int ind_y, const int ind_z) {
		for (int i = 0; i < grid[ind_x][ind_y][ind_z].size(); ++i) {
			if (grid[ind_x][ind_y][ind_z][i] == id) {
				return i;
			}
		}
	}

	long double calculate_h(const int grid_size, const long double h, const std::array<long double, 6> center_pos) {
		const long double element = *(std::max_element(center_pos.begin() + 3, center_pos.end()));
		const long double adj_value = (element * 2) / grid_size;
		return std::fmaxl(adj_value, 1.5 * h);
	}
};

#endif