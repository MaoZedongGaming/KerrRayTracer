#include <vector>

//Vector4d is too AoS for raytracing, we need SoA for better cache performance

// pipeline is relativistic camera tetrad -> shoot photons in camera frame according to the camera's orientation where E = 1 -> project photon local coordinates to global Boyer-Lindquist coordinates -> RKDP integrate first order carter equations -> check if photon escapes, is captured, or hits accretion disk
// escape = sample skybox pixel, capture = photon crosses event horizon = black pixel, accretion disk hit = photon crosses equatorial plane at r > r_isco = sample accretion disk texture pixel + doppler shift + gravitational redshift

enum class PhotonState : uint8_t {
	Active = 0,
	Escaped,
	Captured,
	AccretionDiskHit
};

struct Photons {
	// variable quantities that are continually updated during the simulation
	std::vector <double> t;
	std::vector <double> r;
	std::vector <double> theta;
	std::vector <double> phi;

	// read only impact constant that are set by the camera and the initial conditions of the photon
	std::vector <double> xi;
	std::vector <double> eta;

	// states that are also updated before every integration step, but are not used in the integration itself, branches are bad for performance
	std::vector <float> sign_r; // sign of the radial momentum, +1 for outgoing, -1 for ingoing
	std::vector <float> sign_theta; // sign of the polar momentum, +1 for increasing theta, -1 for decreasing theta

	std::vector <PhotonState> state;
	std::vector <uint32_t> activeIndices; 

	size_t count;

	void reserve(size_t i);
};