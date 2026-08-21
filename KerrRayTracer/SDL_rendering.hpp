#pragma once
#include "parameters.hpp"
#include "general_relativity.hpp"
#include "relativistic_camera.hpp"
#include "SDL3/SDL_pixels.h"
#include "SDL3_image/SDL_image.h"
#include <algorithm>
#include <iostream>
#include <vector>
#include <cmath>
#include <omp.h>

constexpr double MAX_TEMP = 10000.0; //kelvin

constexpr double PI = 3.14159265358979323846;
constexpr double TWO_PI = 6.28318530717958647692;
size_t SKY_WIDTH = 0;
size_t SKY_HEIGHT = 0;

std::vector<uint32_t> unpackSkyField(size_t& outWidth, size_t& outHeight) {
    SDL_Surface* rawSurface = IMG_Load("resources/8k_sun.jpg");
    if (!rawSurface) {
        std::cerr << "couldn't load image: " << SDL_GetError() << "\n";
    }
    SDL_Surface* convertedSurface = SDL_ConvertSurface(rawSurface, SDL_PIXELFORMAT_RGBA32);
    SDL_DestroySurface(rawSurface);

    if (!convertedSurface) {
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

std::vector<uint32_t> skyPixels = unpackSkyField(SKY_WIDTH, SKY_HEIGHT);

uint32_t packRGBA32(uint8_t r, uint8_t g, uint8_t b) {
    return ((uint32_t)r << 24) | ((uint32_t)g << 16) | ((uint32_t)b << 8) | ((uint32_t)255);
}

uint32_t kelvinToRGB(double kelvin) {  // RGBA32 format, Tanner Helland algorithm
    if (kelvin < 1000.0)
        return packRGBA32(0, 0, 0);

    double t = kelvin / 100.0;
    double r, g, b;

    // Calculate Red
    if (t <= 66.0) {
        r = 255.0;
    }
    else {
        r = t - 60.0;
        r = 329.698727446 / std::pow(r, 0.1332047592);
        r = std::max(0.0, std::min(255.0, r));
    }

    // Calculate Green
    if (t <= 66.0) {
        g = t;
        g = 99.4708025861 * std::log(g) - 161.1195634253;
    }
    else {
        g = t - 60.0;
        g = 288.1221695283 / std::pow(g, 0.0755148492);
    }
    g = std::max(0.0, std::min(255.0, g));

    // Calculate Blue
    if (t >= 66.0) {
        b = 255.0;
    }
    else {
        if (t <= 19.0) {
            b = 0.0;
        }
        else {
            b = t - 10.0;
            b = 138.5177312231 * std::log(b) - 305.0447927307;
            b = std::max(0.0, std::min(255.0, b));
        }
    }

    return packRGBA32((uint8_t)r, (uint8_t)g, (uint8_t)b);
}

std::vector<uint32_t> populateTemperature() {
    std::vector<uint32_t> res(40000);
    for (size_t i = 999; i < 40000; ++i)
        res[i] = kelvinToRGB((double)(i));
    return res;
}

std::vector<uint32_t> temperatureLookup = populateTemperature();

uint32_t temperatureToRGB(double t) {
    return temperatureLookup[(size_t)std::clamp(t, 999.0, 40000.0)];
}

double keplerianAngularVelocity(double r) {
	return 1.0 / (r * sqrt(r) + a);
}

double g_factor(double r, double xi) {
    return (sqrt(1.0 - 2 / r + 2 * a / (r * sqrt(r)))) / (1.0 + keplerianAngularVelocity(r) * xi);
}

double novikovThorneTemperature(double r) {
	return MAX_TEMP * pow(r_ISCO / r, 0.75) * pow((1.0 - sqrt(r_ISCO / r)), 0.25);  // Shakura-Sunyaev approximation
}

double observedTemperature(double r, double xi) {
	return g_factor(r, xi) * novikovThorneTemperature(r);
}

uint32_t relativisticBeaming(double r, double xi, uint32_t colour) {
    double intensity = pow(g_factor(r, xi), 3.0);
    uint8_t r0 = (colour >> 24) & 0xFF;
    uint8_t g = (colour >> 16) & 0xFF;
    uint8_t b = (colour >> 8) & 0xFF;
    r0 = (uint8_t)std::clamp(r0 * intensity, 0.0, 255.0);
    g = (uint8_t)std::clamp(g * intensity, 0.0, 255.0);
    b = (uint8_t)std::clamp(b * intensity, 0.0, 255.0);
    return packRGBA32(r0, g, b);
}

uint32_t getAccretionColour(double r, double xi) {
    uint32_t colour = 0;
    colour = temperatureToRGB(observedTemperature(r, xi));
    colour = relativisticBeaming(r, xi, colour);
    return colour;
}

uint32_t sampleSkyField(double theta, double phi) {
    double u = std::fmod(phi, TWO_PI);
    if (u < 0.0)
        u += TWO_PI;
    u /= TWO_PI;
    double v = std::clamp(theta / PI, 0.0, 1.0);
    size_t x = std::clamp((size_t)(u * SKY_WIDTH), 0ull, SKY_WIDTH - 1);
    size_t y = std::clamp((size_t)(v * SKY_HEIGHT), 0ull, SKY_HEIGHT - 1);
    return skyPixels[x + SKY_WIDTH * y];
}

void drawScreen(RelativisticCamera& camera, SDL_Texture* streamTexture) {
    void* pixels = nullptr;
    int pitch = 0;
    #pragma omp parallel for schedule(dynamic, 16)
    for (int photonIndex = 0; photonIndex < camera.width * camera.height; ++photonIndex) {
        switch (camera.photons.state[photonIndex]) {
        case PhotonState::Captured:
        case PhotonState::Active:
            camera.pixelBuffer[photonIndex] = 0;
            break;
        case PhotonState::Escaped:
            camera.pixelBuffer[photonIndex] = sampleSkyField(camera.photons.theta[photonIndex], camera.photons.phi[photonIndex]);
            break;
        case PhotonState::AccretionDiskHit:
            camera.pixelBuffer[photonIndex] = getAccretionColour(camera.photons.r[photonIndex], camera.photons.xi[photonIndex]);
            break;
        }
    }

    SDL_UpdateTexture(streamTexture, nullptr, camera.pixelBuffer.data(), (int)(camera.width * sizeof(uint32_t)));
}