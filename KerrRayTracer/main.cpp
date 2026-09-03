#include <SDL3/SDL.h>
#include <SDL3_Image/SDL_image.h>
#include "SDL_rendering.hpp"
#include <iostream>
#include <omp.h>

int main() {

    #pragma omp parallel
	{
    #pragma omp single
		std::cout << "Running with " << omp_get_num_threads() << " thread(s)\n";
	}
	constexpr int WIDTH = 1920;
	constexpr int HEIGHT = 1080;

	SDL_Init(SDL_INIT_VIDEO);

	SDL_Window* window = SDL_CreateWindow("Kerr Ray Tracer", WIDTH, HEIGHT, 0);
	SDL_Renderer* renderer = SDL_CreateRenderer(window, nullptr);
	
	SDL_Texture* streamTexture = SDL_CreateTexture(renderer, SDL_PIXELFORMAT_RGBA32, SDL_TEXTUREACCESS_TARGET, WIDTH, HEIGHT);

	if (!streamTexture) {
		std::cerr << "failed to create streaming texture, whatcha gunna do: " << SDL_GetError() << "\n";
		return -1;
	}

	bool running = true;

	uint64_t ms_elapsed = 0;
	
	RelativisticCamera camera(WIDTH, HEIGHT);
	camera.setPosition(30.0, 1.5);
	//camera.turnLeft(-0.1);

	camera.initTetrad();
	SDL_Log("generating photons...\n");
	ms_elapsed = SDL_GetTicks();

	camera.generatePhotons();
	ms_elapsed = SDL_GetTicks();
	SDL_Log("%u photons generated in %u ms!\nstarting ray tracing...\n", camera.photons.count, ms_elapsed);

	camera.photons.traceAllRays();
	ms_elapsed = SDL_GetTicks();
	SDL_Log("finished ray tracing in %u ms!\nDrawing pixels to buffer...", ms_elapsed);

	drawScreen(camera, streamTexture);
	ms_elapsed = SDL_GetTicks();
	SDL_Log("finished drawing screen in %u ms!\n", ms_elapsed);

	saveTextureToPNG(renderer, streamTexture, "raymarched_kerr_1080p_no_gamma.png");
	/*
	SDL_Event event;
	while (running) {
		while (SDL_PollEvent(&event)) {
			if (event.type == SDL_EVENT_QUIT)
				running = false;
		}

		camera.initTetrad();
		SDL_Log("generating photons...\n");
		ms_elapsed = SDL_GetTicks();

		camera.generatePhotons();
		ms_elapsed = SDL_GetTicks();
		SDL_Log("%u photons generated in %u ms!\nstarting ray tracing...\n", camera.photons.count, ms_elapsed);

		camera.photons.traceAllRays();
		ms_elapsed = SDL_GetTicks();
		SDL_Log("finished ray tracing in %u ms!\nDrawing pixels to buffer...", ms_elapsed);

		drawScreen(camera, streamTexture);
		ms_elapsed = SDL_GetTicks();
		SDL_Log("finished drawing screen in %u ms!\n", ms_elapsed);

		SDL_RenderClear(renderer);

		SDL_RenderTexture(renderer, streamTexture, nullptr, nullptr);
		SDL_RenderPresent(renderer);
	}*/

	
	SDL_DestroyTexture(streamTexture);
	SDL_DestroyRenderer(renderer);
	SDL_DestroyWindow(window);
	SDL_Quit();
}