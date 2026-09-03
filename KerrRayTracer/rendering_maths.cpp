#pragma once
#include "rendering_maths.hpp"
#include <cstdint>
#include <cmath>
#include <algorithm>

float3& float3::operator += (float3 const& other) {
	r += other.r;
	g += other.g;
	b += other.b;
	return *this;
}

float3& float3::operator *= (float scalar) {
	r *= scalar;
	g *= scalar;
	b *= scalar;
	return *this;
}
float3 float3::operator * (float scalar) const {
	float3 res = *this;
	res *= scalar;
	return res;
}
float3 float3::operator + (float3 const& other) const {
	float3 res = *this;
	res += other;
	return res;
}

float3 operator * (float scalar, float3 const& vec) {
	return vec * scalar;
}

namespace std {
	float3 clamp(float3 const& value, float min, float max) {
		return float3{
			std::clamp(value.r, min, max),
			std::clamp(value.g, min, max),
			std::clamp(value.b, min, max)
		};
	}
}

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

float3 rgbaToFloat3(uint32_t colour) {
	return float3{ getRed(colour) / 255.0f, getGreen(colour) / 255.0f, getBlue(colour) / 255.0f };
}

float3 applyGammaCorrection(float3 colour, float gamma) {
	return float3{ std::pow(colour.r, 1.0f / gamma), std::pow(colour.g, 1.0f / gamma), std::pow(colour.b, 1.0f / gamma) };
}

float3 applyReinhard(float3 colour) {
	return colour * (1.0f / (1.0f + std::max({ colour.r, colour.g, colour.b })));
}

float3 applyACES(float3 colour) {
	float r = colour.r;
	float g = colour.g;
	float b = colour.b;
	r = r * (2.51f * r + 0.03f) / (r * (2.43f * r + 0.59f) + 0.14f);
	g = g * (2.51f * g + 0.03f) / (g * (2.43f * g + 0.59f) + 0.14f);
	b = b * (2.51f * b + 0.03f) / (b * (2.43f * b + 0.59f) + 0.14f);
	return float3{ std::clamp(r, 0.0f, 1.0f), std::clamp(g, 0.0f, 1.0f), std::clamp(b, 0.0f, 1.0f) };
}

uint32_t float3ToRGBA(float3 colour) {
	colour = applyACES(colour);
	//colour = applyGammaCorrection(colour);
	return packRGBA32(
		static_cast<uint8_t>(std::clamp(colour.r * 255.0f, 0.0f, 255.0f)),
		static_cast<uint8_t>(std::clamp(colour.g * 255.0f, 0.0f, 255.0f)),
		static_cast<uint8_t>(std::clamp(colour.b * 255.0f, 0.0f, 255.0f))
	);
}

