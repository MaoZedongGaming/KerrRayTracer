#include <SDL3/SDL.h>
#include <SDL3_Image/SDL_image.h>
#include "SDL_rendering.hpp"
#include <iostream>
#include <omp.h>

bool saveTextureToPNG(SDL_Renderer* renderer, SDL_Texture* texture, const char* filename) {
	// 1. Save current render target to restore it later
	SDL_Texture* previous_target = SDL_GetRenderTarget(renderer);

	// 2. Set target to the texture we want to read
	if (!SDL_SetRenderTarget(renderer, texture)) {
		return false;
	}

	// 3. Get texture dimensions
	float w, h;
	if (!SDL_GetTextureSize(texture, &w, &h)) {
		SDL_SetRenderTarget(renderer, previous_target);
		return false;
	}

	SDL_Rect rect = { 0, 0, (int)w, (int)h };

	// 4. Read pixels into a new SDL3 surface (SDL_RenderReadPixels returns an SDL_Surface* in SDL3)
	SDL_Surface* surface = SDL_RenderReadPixels(renderer, &rect);
	if (!surface) {
		SDL_SetRenderTarget(renderer, previous_target);
		return false;
	}

	// 5. Save image to file (using SDL_SaveBMP or IMG_SavePNG from SDL_image)
	bool success = IMG_SavePNG(surface, filename); // Or SDL_SaveBMP(surface, filename);

	// 6. Clean up
	SDL_DestroySurface(surface);
	SDL_SetRenderTarget(renderer, previous_target);

	return success;
}

int main() {

    #pragma omp parallel
	{
    #pragma omp single
		std::cout << "Running with " << omp_get_num_threads() << " thread(s)\n";
	}
	constexpr int WIDTH = 1000;
	constexpr int HEIGHT = 750;

	SDL_Init(SDL_INIT_VIDEO);

	SDL_Window* window = SDL_CreateWindow("Kerr Ray Tracer", WIDTH, HEIGHT, 0);
	SDL_Renderer* renderer = SDL_CreateRenderer(window, nullptr);
	
	SDL_Texture* streamTexture = SDL_CreateTexture(renderer, SDL_PIXELFORMAT_RGBA32, SDL_TEXTUREACCESS_TARGET, WIDTH, HEIGHT);

	
	//camera.pitch = 0.5;
	if (!streamTexture) {
		std::cerr << "failed to create streaming texture, whatcha gunna do: " << SDL_GetError() << "\n";
		return -1;
	}

	bool running = true;

	uint64_t ms_elapsed = 0;
	
	RelativisticCamera camera(WIDTH, HEIGHT);
	camera.setPosition(30.0, 1.5);

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

	saveTextureToPNG(renderer, streamTexture, "README_2.png");
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