#include "functions.h"
#include "particle.h"
#include <SFML/Graphics.hpp>
#include <SFML/Window.hpp>
#include <iostream>
#include <chrono>

int main()
{
	int i{};
	constexpr int n{ 2 };
	constexpr int pixel_num = 512;
	constexpr float delta_t = .1f;
	short tick{};
	const int compute_ratio = 1;

	//Create the window and particles.
	sf::RenderWindow window(sf::VideoMode(pixel_num, pixel_num), "Particle Simulation");
	sf::View view = window.getDefaultView();
	view.setSize(pixel_num, -pixel_num);
	window.setView(view);
	//window.setFramerateLimit(60);
	//window.setVerticalSyncEnabled(true);

	Particle** part_array = new Particle * [n];
	for (i = 0; i < n; ++i) {
		part_array[i] = new Particle(i, n);
	}
	update_accel(part_array, n, tick);

	//Create shapes via array pointer.
	sf::CircleShape** vertex_array = new sf::CircleShape * [n];
	for (i = 0; i < n; ++i) {
		vertex_array[i] = new sf::CircleShape[n];
		(*vertex_array[i]).setRadius(2.f);
		(*vertex_array[i]).setPointCount(6);
	}

	/*float alpha{-1.218f}, beta{.783f}, gamma{.209f};*/
	float alpha{}, beta{}, gamma{}, scale{};
	const float x_cam{}, y_cam{};
	float display_radius{ 1e8 };

	//Particle calculations with velocity Verlet and timer.
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
					display_radius -= 5e6;
					if (display_radius < 0)
						display_radius = 0;
				}
				if (event.key.code == sf::Keyboard::K) {
					display_radius += 5e6;
				}
			}
		}

		//Start of main window loop.
		window.clear(sf::Color::Black);

		update_pos(part_array, n, delta_t);
		update_accel(part_array, n, tick);
		update_vel(part_array, n, delta_t);

		if (tick % compute_ratio == 0) {
			//print_Energy(part_array, n);

			//Rendering functions.
			update_vertex_pos(part_array, vertex_array, n, pixel_num, display_radius, alpha, beta, gamma, tick);
			for (i = 0; i < n; ++i) {
				window.draw((*vertex_array[i]));
			}
			window.display();
		}

		//if (tick == SHRT_MAX) {
		//	auto t2 = std::chrono::high_resolution_clock::now();
		//	auto ms_int = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
		//	std::cout << ms_int.count() << '\n';
		//}

		tick += 1;

	}
	delete[] part_array;
	delete[] vertex_array;
}

