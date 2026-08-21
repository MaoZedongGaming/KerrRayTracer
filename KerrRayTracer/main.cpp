#include <SDL3/SDL.h>
#include <SDL3_Image/SDL_image.h>
#include "SDL_rendering.hpp"
#include <iostream>

int main() {
	constexpr int WIDTH = 800;
	constexpr int HEIGHT = 600;

	SDL_Init(SDL_INIT_VIDEO);

	SDL_Window* window = SDL_CreateWindow("Kerr Ray Tracer", WIDTH, HEIGHT, 0);
	SDL_Renderer* renderer = SDL_CreateRenderer(window, nullptr);
	
	SDL_Texture* streamTexture = SDL_CreateTexture(renderer, SDL_PIXELFORMAT_RGBA32, SDL_TEXTUREACCESS_STREAMING, WIDTH, HEIGHT);

	RelativisticCamera camera(WIDTH, HEIGHT);
	camera.setPosition(50.0, 1.0);

	if (!streamTexture) {
		std::cerr << "failed to create streaming texture, whatcha gunna do: " << SDL_GetError() << "\n";
		return -1;
	}

	bool running = true;
	SDL_Event event;
	while (running) {
		while (SDL_PollEvent(&event)) {
			if (event.type == SDL_EVENT_QUIT)
				running = false;
		}
		
		camera.initTetrad();
		SDL_Log("generating photons...\n");
		camera.generatePhotons();
		SDL_Log("%llu photons generated!\nstarting ray tracing...\n", camera.photons.count);
		camera.photons.traceAllRays();
		SDL_Log("finished ray tracing!\n");
		drawScreen(camera, streamTexture);
		//SDL_RenderClear(renderer);
		SDL_FRect dstRect = { 0.0f, 0.0f, static_cast<float>(WIDTH), static_cast<float>(HEIGHT) };
		SDL_RenderTexture(renderer, streamTexture, nullptr, &dstRect);
		SDL_RenderPresent(renderer);
	}

	
	SDL_DestroyTexture(streamTexture);
	SDL_DestroyRenderer(renderer);
	SDL_DestroyWindow(window);
	SDL_Quit();
}