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

	// constant quantities that are set by the camera and the initial conditions of the photon
	std::vector <double> E; // E = - p_t, the energy of the photon
	std::vector <double> L_z; // L_z = p_phi, the angular momentum of the photon about the axis of rotation
	std::vector <double> Q; // Carter constant, a conserved quantity related to the total angular momentum of the photon

	// states that are also updated before every integration step, but are not used in the integration itself, branches are bad for performance
	std::vector <double> sign_r; // sign of the radial momentum, +1 for outgoing, -1 for ingoing
	std::vector <double> sign_theta; // sign of the polar momentum, +1 for increasing theta, -1 for decreasing theta

	std::vector <PhotonState> state; // the state of the photon, whether it is active, escaped, captured, or hit the accretion disk
	std::vector <uint32_t> activeIndices; // the indices of the active photons, used for efficient integration of only the active photons

	size_t count;
};

// float precision is enough for raytracing, we don't need double precision, and float is faster and uses less memory, which is important for large numbers of photons
// nah double precision for photons, float precision for camera and accretion disk, because photons are the most important part of the simulation and we want to avoid numerical errors in their trajectorie