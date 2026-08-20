
#include "maths.hpp"
#include "photons.hpp"

// the camera object has to be a relativistic object within the kerr spacetime, no god views allowed compared to schwarzschild spacetime, and for the E = 1 generalisation to work the frame must be stationary so ZAMO frame
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
	Photons photons;

	float fov = 60.0f; // field of view in degrees
	int width;
	int height;

	double pitch = 0.0f; // camera pitch in degrees
	double yaw = 0.0f; // camera yaw in degrees

	void initTetrad();
	void setPosition(double r, double theta);
	void generatePhotons();

	RelativisticCamera(int w, int h);
};