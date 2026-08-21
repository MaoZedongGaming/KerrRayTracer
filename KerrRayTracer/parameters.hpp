#pragma once 
#include <cmath>

constexpr double a = 0.6;

constexpr double r_sky = 200.0; // skybox radius from the centre

constexpr double r_acceretion = 20.0; // maximal radius of acceretion disk

namespace {
	double calculateOuterHorizon() {
		return 1.0 + sqrt(1 - a * a);
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