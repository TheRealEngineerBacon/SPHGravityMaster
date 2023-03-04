#ifndef SMOOTHGRID_H
#define SMOOTHGRID_H
#include "particle.h"
#include <vector>
#include <array>


class SmoothGrid {
public:
	std::vector<std::vector<std::vector<std::vector<int>>>> grid;

	SmoothGrid(int size, long double h) {
		grid.resize(size, std::vector<std::vector<std::vector<int>>>
				   (size, std::vector<std::vector<int>>
				   (size, std::vector<int>())));
	}

	void on_grid(Particle** array, int n, int grid_size, long double h, std::array<long double, 6> center_pos) {
		long double h_adj = calculate_h(grid_size, h, center_pos);
		for (int i = 0; i < n; ++i) {
			if ((*array[i]).on_grid == true) {
				int offset = floor(grid_size / 2.f);
				int a{}, b{}, c{};
				long double x{ (*array[i]).x }, y{ (*array[i]).y }, z{ (*array[i]).z };
				long double h_half = h_adj / 2;
				a = floor((x + h_half - center_pos[0]) / h_adj) + offset;
				b = floor((y + h_half - center_pos[1]) / h_adj) + offset;
				c = floor((z + h_half - center_pos[2]) / h_adj) + offset;
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

	void place(Particle** array, int i, int grid_size, long double h, std::array<long double, 6> center_pos, long double h_adj) {
		int offset = floor(grid_size / 2.f);
		int a{}, b{}, c{};
		long double x{ (*array[i]).x }, y{ (*array[i]).y }, z{ (*array[i]).z };
		long double h_half = h_adj / 2;
		a = floor((x + h_half - center_pos[0]) / h_adj) + offset;
		b = floor((y + h_half - center_pos[1]) / h_adj) + offset;
		c = floor((z + h_half - center_pos[2]) / h_adj) + offset;
		grid[a][b][c].push_back((*array[i]).id);
		(*array[i]).ind_x = a;
		(*array[i]).ind_y = b;
		(*array[i]).ind_z = c;
	}

	int find_element(int id, int ind_x, int ind_y, int ind_z) {
		for (int i = 0; i < grid[ind_x][ind_y][ind_z].size(); ++i) {
			if (grid[ind_x][ind_y][ind_z][i] == id) {
				return i;
			}
		}
	}

	long double calculate_h(int grid_size, long double h, std::array<long double, 6> center_pos) {
		long double max_element = *std::max_element(center_pos.begin() + 3, center_pos.end());
		int adj_value = (max_element * 2) / grid_size;
		return std::fmaxl(adj_value, h);
	}

	std::vector<int> return_surr(Particle** array, int focus, int grid_size) {
		std::vector<int> temp;
		int sindex_x{ (*array[focus]).ind_x }, sindex_y{ (*array[focus]).ind_y }, sindex_z{ (*array[focus]).ind_z };
		for (int i = sindex_x - 1; i < (sindex_x + 2); ++i) {
			if (i < 0 || i > grid_size - 1)
				continue;
			for (int j = sindex_y - 1; j < (sindex_y + 2); ++j) {
				if (j < 0 || j > grid_size - 1)
					continue;
				for (int k = sindex_z - 1; k < (sindex_z + 2); ++k) {
					if (k < 0 || k > grid_size - 1)
						continue;
					for (int p = 0; p < grid[i][j][k].size(); ++p)
						temp.push_back(grid[i][j][k][p]);
				}
			}
		}
		return temp;
	}
};


#endif