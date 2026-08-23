#pragma once
#include <vector>

//Vector4d is too AoS for raytracing, we need SoA for better cache performance

// pipeline is relativistic camera tetrad -> shoot photons in camera frame according to the camera's orientation where E = 1 -> project photon local coordinates to global Boyer-Lindquist coordinates -> RKDP integrate first order carter equations -> check if photon escapes, is captured, or hits accretion disk
// escape = sample skybox pixel, capture = photon crosses event horizon = black pixel, accretion disk hit = photon crosses equatorial plane at r > r_isco = sample accretion disk texture pixel + doppler shift + gravitational redshift

enum class PhotonState : uint8_t {
	Active = 0,
	Escaped,
	Captured,
	AccretionDiskHit,
	PhotonSphere
};

struct PhotonDerivative {
	//double dt;
	double dr;
	double dtheta;
	double dphi;
};

PhotonDerivative evaluate(double r, double theta, float sign_r, float sign_theta, double xi, double eta);

struct Photons {
	// variable quantities that are continually updated during the simulation
	//std::vector <double> t;  // I don't think you'd ever need coordinate time for anything in rendering 
	std::vector <double> r;
	std::vector <double> theta;
	std::vector <double> phi;

	// you usually don't need to store or cache the momenta but this is the k1 necessaary for RKDP's FSAL, new k1 = old k7, also they must be evaluated as carter's equations and not the initial condition direction vector
	std::vector <double> dr;
	std::vector <double> dtheta;
	std::vector <double> dphi;

	// read only impact constant that are set by the camera and the initial conditions of the photon
	std::vector <double> xi;  // L_Z / E
	std::vector <double> eta;  // Q / E^2

	// states that are also updated before every integration step, but are not used in the integration itself, branches are bad for performance
	std::vector <float> sign_r; // sign of the radial momentum, +1 for outgoing, -1 for ingoing
	std::vector <float> sign_theta; // sign of the polar momentum, +1 for increasing theta, -1 for decreasing theta

	std::vector <PhotonState> state;
	//std::vector <size_t> activeIndices; 

	std::vector <double> dlambda;  // per photon adaptative time step
	// would using add a k1 array for FSAL in RKDP but that would be bad for the GPU

	size_t count = 0;

	//void reserve(size_t i);
	void resize(size_t i);
	//void clear();
	// adaptive RKDP, mino time, impact constants, affine parameter starts at 0.01
	void rkdpStepRay(size_t i);
	void traceAllRays();
};