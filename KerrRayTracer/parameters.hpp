#pragma once 
#include <cmath>

constexpr double PI = 3.14159265358979323846;
constexpr double TWO_PI = 6.28318530717958647692;
constexpr double PI_2 = 1.57079632679;

template <typename T>
constexpr int sgn(T x) {
	return ((x > T(0)) - (x < T(0)));
}

template <>
constexpr int sgn(double x) {
	return ((x > 1e-14) - (1e-14 > x));
}

constexpr double a = 0.6;

constexpr double r_sky = 40.0; // skybox radius from the centre

constexpr double r_acceretion = 20.0; // maximal radius of acceretion disk


constexpr double MAX_TEMP = 4000.0; //kelvin


namespace {
	double calculateOuterHorizon() {
		return 1.0 + sqrt(1.0 - a * a);
	}
	double Z1() {
		return 1.0 + cbrt(1.0 - a * a) * (cbrt(1.0 + a) + cbrt(1.0 - a));
	}
	double Z2() {
		return sqrt(3 * a * a + Z1() * Z1());
	}
	double calculateCorotatingISCORadius() {
		return 3.0 + Z2() - sqrt((3.0 - Z1()) * (3.0 + Z1() + 2.0 * Z2()));
	}
}

const double r_horizon = calculateOuterHorizon();  // outer event horizon radius, no constexpr sqrt before C++26 :(

const double r_ISCO = calculateCorotatingISCORadius(); // radius of the innermost stable circular orbit for corotating/prograde gasses 