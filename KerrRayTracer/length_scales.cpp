#include "length_scales.hpp"
#include <cmath>

constexpr double sigma(double radius, double theta) {
	return radius * radius + a * a * cos(theta) * cos(theta);
}

constexpr double delta(double radius) {
	return radius * radius - 2 * radius + a * a;
}

// wtf am I doing