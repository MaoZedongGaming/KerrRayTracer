#pragma once
#include <cstdint>

struct float3 {
	float r = 0.0f;
	float g = 0.0f;
	float b = 0.0f;
	float3& operator += (float3 const& other);
	float3& operator *= (float scalar);
	float3 operator * (float scalar) const;
	float3 operator + (float3 const& other) const;
};

float3 operator * (float scalar, float3 const& vec);

namespace std {
	float3 clamp(float3 const& value, float min, float max);
}

constexpr uint32_t packRGBA32(uint8_t r, uint8_t g, uint8_t b);

constexpr uint8_t getRed(uint32_t colour);

constexpr uint8_t getGreen(uint32_t colour);

constexpr uint8_t getBlue(uint32_t colour);

float3 rgbaToFloat3(uint32_t colour);

float3 applyGammaCorrection(float3 colour, float gamma = 2.2f);

float3 applyReinhard(float3 colour);

float3 applyACES(float3 colour);

uint32_t float3ToRGBA(float3 colour);
