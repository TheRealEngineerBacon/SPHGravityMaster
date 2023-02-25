#include "constants.h"
#include "functions.h"
#include "particle.h"
#include <SFML/Graphics.hpp>
#include <SFML/Window.hpp>
#include <iostream>
#include <chrono>

int main()
{
	//Create the window and particles.
	sf::RenderWindow window(sf::VideoMode(pixel_num, pixel_num), "Particle Simulation");
	sf::View view = window.getDefaultView();
	view.setSize(pixel_num, -pixel_num);
	window.setView(view);
	//window.setFramerateLimit(60);
	//window.setVerticalSyncEnabled(true);

	Particle** part_array = new Particle * [n];
	for (int i = 0; i < n; ++i) {
		part_array[i] = new Particle(i, n);
	}
	update_grav(part_array, n, tick);

	//Create shapes via array pointer.
	sf::CircleShape** vertex_array = new sf::CircleShape * [n];
	for (int i = 0; i < n; ++i) {
		vertex_array[i] = new sf::CircleShape[n];
		(*vertex_array[i]).setRadius(2.f);
		(*vertex_array[i]).setPointCount(8);
	}

	//std::chrono::steady_clock::time_point t1, t2;
	auto t1 = std::chrono::high_resolution_clock::now();
	
	//Create main window loop.
	while (window.isOpen()) {
		window.setKeyRepeatEnabled(true);
		sf::Event event;

		//Event handling.
		while (window.pollEvent(event)) {
			switch (event.type) {
			case sf::Event::Closed:
				window.close();
			case sf::Event::KeyPressed:
				if (event.key.code == sf::Keyboard::LShift) {
					std::cout << "Closing key pressed.\n";
					window.close();
				}
				if (event.key.code == sf::Keyboard::R) {
					alpha = 0, beta = 0, gamma = 0;
				}
				if (event.key.code == sf::Keyboard::Left) {
					beta -= 0.087f;
				}
				if (event.key.code == sf::Keyboard::Right) {
					beta += 0.087f;
				}
				if (event.key.code == sf::Keyboard::Down) {
					alpha += 0.087f;
				}
				if (event.key.code == sf::Keyboard::Up) {
					alpha -= 0.087f;
				}
				if (event.key.code == sf::Keyboard::Delete) {
					gamma += 0.087f;
				}
				if (event.key.code == sf::Keyboard::End) {
					gamma -= 0.087f;
				}
				if (event.key.code == sf::Keyboard::I) {
					display_radius -= 1e6;
					if (display_radius < 0)
						display_radius = 0;
				}
				if (event.key.code == sf::Keyboard::K) {
					display_radius += 1e6;
				}
				if (event.key.code == sf::Keyboard::U) {
					focus += 1;
					if (focus > (n - 1)) {
					focus = n - 1;
					}
					std::cout << focus << '\n';
				}
				if (event.key.code == sf::Keyboard::J) {
					focus -= 1;
					if (focus < 0) {
						focus = 0;
					}
					std::cout << focus << '\n';
				}
				if (event.key.code == sf::Keyboard::W) {
					y_mov -= 1e6f;
				}
				if (event.key.code == sf::Keyboard::S) {
					y_mov += 1e6f;
				}
				if (event.key.code == sf::Keyboard::A) {
					x_mov += 1e6f;
				}
				if (event.key.code == sf::Keyboard::D) {
					x_mov -= 1e6f;
				}
				if (event.key.code == sf::Keyboard::Space) {
					if (pause_state == false) {
						delta_t = 0.0;
						pause_state = true;
					}
					else {
						delta_t = 0.5;
						pause_state = false;
					}
				}

			}
		}

		//Start of main window loop.
		window.clear(sf::Color::Black);

		update_pos(part_array, n, delta_t);
		update_grav(part_array, n, tick);
		update_fluid(part_array, n, tick);
		update_vel(part_array, n, delta_t);
		update_temp(part_array, n, delta_t);

		if (tick % compute_ratio == 0) {

			//Rendering functions.
			update_vertex_pos(part_array, vertex_array, n, pixel_num, display_radius, alpha, beta, gamma, tick, focus, x_mov, y_mov);
			for (int i = 0; i < n; ++i) {
				window.draw((*vertex_array[i]));
			}
			window.display();
		}

		if (tick % 8000 == 0) {
			std::cout << std::setprecision(12) << (*part_array[focus]).temp_f << '\n';
		}

		int image_num{};
		if (tick == SHRT_MAX) {
			/*auto t2 = std::chrono::high_resolution_clock::now();
			auto ms_int = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
			std::cout << ms_int.count() << '\n';*/
			//print_Data(part_array, n, focus);
			check_Singularity(part_array, n);
		}

		tick += 1;
	}
	delete[] part_array;
	delete[] vertex_array;
}

