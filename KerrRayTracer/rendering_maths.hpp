#pragma once
#include <array>
#include <random>
#include <cmath>
#include <algorithm>

using float3 = std::array<float, 3>;

constexpr uint32_t packRGBA32(uint8_t r, uint8_t g, uint8_t b) {
    return ((uint32_t)255 << 24) | ((uint32_t)b << 16) | ((uint32_t)g << 8) | ((uint32_t)r);
}

constexpr uint8_t getRed(uint32_t colour) {
    return (colour & 0xFF);
}

constexpr uint8_t getGreen(uint32_t colour) {
    return ((colour >> 8) & 0xFF);
}

constexpr uint8_t getBlue(uint32_t colour) {
    return ((colour >> 16) & 0xFF);
}

constexpr float3 rgbaToFloat3(uint32_t colour) {
    return float3{ getRed(colour) / 255.0f, getGreen(colour) / 255.0f, getBlue(colour) / 255.0f };
}

constexpr uint32_t float3ToRGBA(float3 fcolour) {
    return packRGBA32((uint8_t)std::clamp(fcolour[0] * 255.0f, 0.0f, 255.0f), (uint8_t)std::clamp(fcolour[1] * 255.0f, 0.0f, 255.0f), (uint8_t)std::clamp(fcolour[2] * 255.0f, 0.0f, 255.0f));
}

//double applyReinhard(double x) {
//    return x / (1.0 + x);
//}
//
//double applyACES(double x) {
//    return x * (2.51 * x + 0.03) / (x * (2.43 * x + 0.59) + 0.14);
//}