#include "maths.hpp"

// the camera object has to be a relativistic object within the kerr spacetime, no god views allowed compared to schwarzschild spacetime

struct Photons;

struct Tetrad {
	Vector4d e0; // timelike vector, camera 4-velocity
	Vector4d e1; // spacelike vector, camera right direction
	Vector4d e2; // spacelike vector, camera up direction
	Vector4d e3; // spacelike vector, camera forward direction
};

class RelativisticCamera {
public:
	Tetrad frame; // camera orthonormal tetrad
	Vector4d position; // (t, r, theta, phi) in Boyer-Lindquist coordinates
	Vector4d velocity; // (dt/dtau, dr/dtau, dtheta/dtau, dphi/dtau) in Boyer-Lindquist coordinates

	float fov = 60.0f; // field of view in degrees
	int width;
	int height;

	double roll = 0.0f; // camera roll in degrees
	double pitch = 0.0f; // camera pitch in degrees
	double yaw = 0.0f; // camera yaw in degrees

	void initTetrad();
	void setPosition(Vector4d const& pos);
	void generatePhotons(Photons& photons);
};