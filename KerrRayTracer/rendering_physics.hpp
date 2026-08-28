#pragma once
#include "rendering_maths.hpp"
#include "parameters.hpp"
#include "config.hpp"
#include <vector>
#include <random>
#include <cmath>
#include <algorithm>

// copy pasted I can't bother to implement brownian noise myself
class PrecomputedfBm2D {
public:
    explicit PrecomputedfBm2D(size_t tableSize = 512, unsigned int seed = std::random_device{}())
        : m_mask(tableSize - 1) {
        // Ensure table size is a power of two for lightning-fast bitwise masking
        if ((tableSize & m_mask) != 0) {
            throw std::invalid_argument("Table size must be a power of two.");
        }

        m_noiseTable.resize(tableSize * tableSize);
        std::mt19937 gen(seed);
        std::uniform_real_distribution<float> dist(-1.0f, 1.0f);

        // Populate the look-up table with uniform white noise
        std::generate(m_noiseTable.begin(), m_noiseTable.end(), [&]() { return dist(gen); });
    }

    // High-performance 2D fBm Noise evaluation function
    [[nodiscard]] float sample(float x, float y, int octaves = 4, float lacunarity = 2.0f, float gain = 0.5f) const noexcept {
        float total = 0.0f;
        float amplitude = 1.0f;
        float maxValue = 0.0f; // Track maximum possible value for clean normalization

        for (int i = 0; i < octaves; ++i) {
            total += amplitude * getNoise(x, y);
            maxValue += amplitude;

            x *= lacunarity;
            y *= lacunarity;
            amplitude *= gain;
        }

        return total / maxValue; // Returns normalized value clamped between [-1.0, 1.0]
    }

private:
    // Quintic blending curve (S-curve) for ultra-smooth edge blending
    [[nodiscard]] static constexpr float fade(float t) noexcept {
        return t * t * t * (t * (t * 6.0f - 15.0f) + 10.0f);
    }

    // Linear interpolation
    [[nodiscard]] static constexpr float lerp(float t, float a, float b) noexcept {
        return a + t * (b - a);
    }

    // Core value look-up and interpolation logic
    [[nodiscard]] float getNoise(float x, float y) const noexcept {
        // Find grid cell coordinates using floor
        auto xFloor = static_cast<int>(std::floor(x));
        auto yFloor = static_cast<int>(std::floor(y));

        // Get fractional offsets inside the specific grid cell [0.0, 1.0)
        float xFract = x - static_cast<float>(xFloor);
        float yFract = y - static_cast<float>(yFloor);

        // Apply bitwise masking to tile coordinates smoothly across bounds
        size_t x0 = static_cast<size_t>(xFloor) & m_mask;
        size_t x1 = (x0 + 1) & m_mask;
        size_t y0 = static_cast<size_t>(yFloor) & m_mask;
        size_t y1 = (y0 + 1) & m_mask;

        // Fetch precomputed corner values from 1D flat array mapping
        float c00 = m_noiseTable[y0 * (m_mask + 1) + x0];
        float c10 = m_noiseTable[y0 * (m_mask + 1) + x1];
        float c01 = m_noiseTable[y1 * (m_mask + 1) + x0];
        float c11 = m_noiseTable[y1 * (m_mask + 1) + x1];

        // Smooth out coordinates using fade curve weights
        float u = fade(xFract);
        float v = fade(yFract);

        // Perform bilinear blend across 4 corners
        return lerp(v, lerp(u, c00, c10), lerp(u, c01, c11));
    }

    std::vector<float> m_noiseTable;
    size_t m_mask;
};

namespace {
    // RGBA32 format, Tanner Helland algorithm
    uint32_t kelvinToRGB(double kelvin) {
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
        std::vector<uint32_t> res(40000, packRGBA32(0, 0, 0));
        for (size_t i = 1000; i < 40000; ++i)
            res[i] = kelvinToRGB((double)(i));
        return res;
    }
}

// convert to textures for hlsl!!!

PrecomputedfBm2D noiseGen(512);

std::vector<uint32_t> temperatureLookup = populateTemperature();

uint32_t temperatureToRGB(double t) {
    return temperatureLookup[(size_t)std::clamp(t, 0.0, 40000.0)];
}

double keplerianAngularVelocity(double r) {
    return 1.0 / (r * sqrt(r) + a);
}

double diskRotationAngle(double r, double phi, double cameraTime) {
    return phi - keplerianAngularVelocity(r) * cameraTime;
}

float diskDensity(double r, double phi, double cameraTime) {
    double phi_disk = diskRotationAngle(r, phi, cameraTime);
    return noiseGen.sample((float)(r * cos(phi_disk)), (float)(r * sin(phi_disk)));
}

//bool intersectAccretionDisk(double r, double theta) {
//	return ((abs(r * cos(theta)) <= accretionThickness(r) + 1e-3) && (r_ISCO <= r && r <= r_acceretion));
//}

constexpr bool crossedEquatorialPlane(double theta_old, double theta_new) {
    return (std::min(theta_old, theta_new) <= PI_2 && PI_2 <= std::max(theta_old, theta_new));
}

bool intersectAccretionDisk(double r, double theta_old, double theta_new) {
    return crossedEquatorialPlane(theta_old, theta_new) && (r_ISCO < r && r <= r_acceretion);
}

double g_factor(double r, double xi) {
    return (sqrt(1.0 - 2.0 / r + 2.0 * a / (r * sqrt(r)))) / (1.0 + keplerianAngularVelocity(r) * xi);
}

// Shakura-Sunyaev approximation
double novikovThorneTemperature(double r) {
    return MAX_TEMP * pow((r_ISCO / r), 0.75) * sqrt(sqrt(1.0 - sqrt(r_ISCO / r)));
}

double observedTemperature(double r, double xi) {
    return g_factor(r, xi) * novikovThorneTemperature(r);
}

uint32_t applyDensity(double r, double phi, double cameraTime, uint32_t colour) {
    float density = diskDensity(r, phi, cameraTime);
    uint8_t r0 = getRed(colour);
    uint8_t g = getGreen(colour);
    uint8_t b = getBlue(colour);
    r0 = (uint8_t)std::clamp(r0 * density, 0.0f, 255.0f);
    g = (uint8_t)std::clamp(g * density, 0.0f, 255.0f);
    b = (uint8_t)std::clamp(b * density, 0.0f, 255.0f);
    return packRGBA32(r0, g, b);
}

float3 applyDensity(double r, double phi, double cameraTime, float3 colour) {
    float density = diskDensity(r, phi, cameraTime);
    return float3{ std::clamp(colour[0] * density, 0.0f, 1.0f), std::clamp(colour[1] * density, 0.0f, 1.0f), std::clamp(colour[2] * density, 0.0f, 1.0f) };
}

uint32_t relativisticBeaming(double r, double xi, uint32_t colour) {
    double intensity = pow(g_factor(r, xi), 4.0);
    //double reinhard_intensity = intensity / (1.0 + intensity);  // looks hideous
    //double ACES_intensity = intensity * (2.51 * intensity + 0.03) / (intensity * (2.43 * intensity + 0.59) + 0.14);  // both look pretty bad
    uint8_t r0 = getRed(colour);
    uint8_t g = getGreen(colour);
    uint8_t b = getBlue(colour);
    r0 = (uint8_t)std::clamp(r0 * intensity, 0.0, 255.0);
    g = (uint8_t)std::clamp(g * intensity, 0.0, 255.0);
    b = (uint8_t)std::clamp(b * intensity, 0.0, 255.0);
    return packRGBA32(r0, g, b);
}

float3 relativisticBeaming(double r, double xi, float3 colour) {
    float intensity = (float)pow(g_factor(r, xi), 4.0);
    return float3{ std::clamp(colour[0] * intensity, 0.0f, 1.0f), std::clamp(colour[1] * intensity, 0.0f, 1.0f), std::clamp(colour[2] * intensity, 0.0f, 1.0f) };
}

uint32_t getAccretionColour(double r, double xi, double phi = 0, double cameraTime = 0) {
    uint32_t colour = packRGBA32(0, 0, 0);
    colour = temperatureToRGB(observedTemperature(r, xi));
    if constexpr (ENABLE_DOPPLER_BEAMING) {
        colour = relativisticBeaming(r, xi, colour);
    }
    if constexpr (ENABLE_DISK_NOISE) {
        colour = applyDensity(r, phi, cameraTime, colour);
    }
    return colour;
}