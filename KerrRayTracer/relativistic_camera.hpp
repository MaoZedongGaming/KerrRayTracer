#pragma once
#include "maths.hpp"
#include "photons.hpp"
#include "parameters.hpp"
#include <vector>

// the camera object has to be a relativistic object within the kerr spacetime, no god views allowed compared to schwarzschild spacetime, and for the E = 1 generalisation to work the frame must be stationary so ZAMO frame
struct Tetrad {
	Vector4d e0; // timelike vector, camera 4-velocity
	Vector4d e1; // spacelike vector, camera right direction
	Vector4d e2; // spacelike vector, camera up direction
	Vector4d e3; // spacelike vector, camera forward direction
};

class RelativisticCamera {
public:
	Photons photons;
	Tetrad frame; // camera orthonormal tetrad
	Vector4d position; // (t, r, theta, phi) in Boyer-Lindquist coordinates
	Vector4d velocity; // (dt/dtau, dr/dtau, dtheta/dtau, dphi/dtau) in Boyer-Lindquist coordinates
	std::vector<uint32_t> pixelBuffer;

	float fov = 60.0f * (float)PI / 180.0f; // field of view in radians
	size_t width;
	size_t height;

	double pitch = 0.0; // camera pitch in degrees
	double yaw = 0.0; // camera yaw in degrees

	void initTetrad();
	void setPosition(double r, double theta, double phi = 0.0);
	void generatePhotons();
	void turnLeft(double y);
	void turnDown(double p);

	RelativisticCamera(size_t w, size_t h);
};