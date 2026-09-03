#pragma once
#include "parameters.hpp"
#include "general_relativity.hpp"
#include "relativistic_camera.hpp"
#include "rendering_maths.hpp"
#include "config.hpp"
#include "SDL3/SDL_pixels.h"
#include "SDL3_image/SDL_image.h"
#include <algorithm>
#include <iostream>
#include <vector>
#include <cmath>
#include <omp.h>

size_t SKY_WIDTH = 0;
size_t SKY_HEIGHT = 0;

std::vector<uint32_t> unpackImage(char const* filepath, size_t& outWidth, size_t& outHeight) {
    SDL_Surface* rawSurface = IMG_Load(filepath);
    if (rawSurface == nullptr) {
        std::cerr << "couldn't load image: " << SDL_GetError() << "\n";
    }
    SDL_Surface* convertedSurface = SDL_ConvertSurface(rawSurface, SDL_PIXELFORMAT_RGBA32);
    SDL_DestroySurface(rawSurface);

    if (convertedSurface == nullptr) {
        std::cerr << "couldn't convert surface: " << SDL_GetError() << "\n";
        SDL_Quit();
    }

    outWidth = convertedSurface->w;
    outHeight = convertedSurface->h;
    std::vector<uint32_t> res(outWidth * outHeight);

    uint8_t* sourcePixels = static_cast<uint8_t*>(convertedSurface->pixels);
    int sourcePitch = convertedSurface->pitch;

    for (int y = 0; y < outHeight; ++y) {
        uint32_t* rowStart = reinterpret_cast<uint32_t*>(sourcePixels + (y * sourcePitch));
        uint32_t* destStart = &res[y * outWidth];
        std::memcpy(destStart, rowStart, outWidth * sizeof(uint32_t));
    }
    
    SDL_DestroySurface(convertedSurface);
    return res;
}

// LUTs, convert to textures when porting to hlsl

std::vector<uint32_t> skyPixels = unpackImage("resources/2k_stars.png", SKY_WIDTH, SKY_HEIGHT);

uint32_t sampleSkyField(double theta, double phi) {
    double u = phi / TWO_PI;
    double v = std::clamp(theta / PI, 0.0, 1.0);
    size_t x = (size_t)(u * SKY_WIDTH) % SKY_WIDTH;
    size_t y = std::clamp((size_t)(v * SKY_HEIGHT), 0ull, SKY_HEIGHT - 1);
    return skyPixels[x + SKY_WIDTH * y];
}

void drawScreen(RelativisticCamera& camera, SDL_Texture* streamTexture) {
    #pragma omp parallel for schedule(dynamic, 16)
    for (int i = 0; i < camera.width * camera.height; ++i) {
        switch (camera.photons.state[i]) {
        case PhotonState::Active:
            if constexpr (PIXEL_DEBUG) {
                 camera.pixelBuffer[i] = packRGBA32(255, 0, 255);
                 break;
            }
			[[fallthrough]];
        case PhotonState::Captured:
            camera.pixelBuffer[i] = float3ToRGBA(camera.photons.accumulatedColour[i]);
            break;
        case PhotonState::Escaped:
            if constexpr (ENABLE_OPAQUE_DISK) {
                camera.pixelBuffer[i] = sampleSkyField(camera.photons.theta[i], camera.photons.phi[i]);
            }
            if constexpr (!ENABLE_OPAQUE_DISK) {
                float3 skyboxColour = rgbaToFloat3(sampleSkyField(camera.photons.theta[i], camera.photons.phi[i]));
				float3 gasColour = camera.photons.accumulatedColour[i] + camera.photons.transmittance[i] * skyboxColour * 0.7f;  // extra factor to dim the skyfield to emphasise the accretion disk
                camera.pixelBuffer[i] = float3ToRGBA(gasColour);
            }
            break;
        case PhotonState::AccretionDiskHit:
            camera.pixelBuffer[i] = float3ToRGBA(camera.photons.accumulatedColour[i]);
            break;
        }
    }

    SDL_UpdateTexture(streamTexture, nullptr, camera.pixelBuffer.data(), (int)(camera.width * sizeof(uint32_t)));
}

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