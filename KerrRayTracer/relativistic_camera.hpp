#include "maths.hpp"

struct Tetrad {
	Vector4d e0; // timelike vector, camera 4-velocity
	Vector4d e1; // spacelike vector, camera right direction
	Vector4d e2; // spacelike vector, camera up direction
	Vector4d e3; // spacelike vector, camera forward direction
};

class RelativisticCamera {
public:
	Vector4d position; // (t, r, theta, phi) in Boyer-Lindquist coordinates
	Vector4d velocity; // (dt/dtau, dr/dtau, dtheta/dtau, dphi/dtau) in Boyer-Lindquist coordinates
	Tetrad frame; // camera orthonormal tetrad

	float fov = 60.0f; // field of view in degrees
	int width;
	int height;
};