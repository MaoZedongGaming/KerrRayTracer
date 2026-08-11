#include "length_scales.hpp"
#include <cmath>

constexpr double spinParameter(double spin, double mass) {
	return spin / (mass * mass);
}

constexpr double schwarzschildRadius(double mass) {
	return 2 * mass;
}

constexpr double sigma(double radius, double theta, double mass, double spin) {
	return radius * radius + spinParameter(spin, mass) * spinParameter(spin, mass) * cos(theta) * cos(theta);
}

constexpr double delta(double radius, double mass, double spin) {
	return radius * radius - 2 * mass * radius + spinParameter(spin, mass) * spinParameter(spin, mass);
}