#include "maths.hpp"
#include <vector>

struct Photons {
	std::vector<Vector<double, 4>> positions; // (t, r, theta, phi)
	std::vector<Vector<double, 4>> momenta; // (E, p_r, p_theta, p_phi)
	// E = -p_t
	// p_r = dr/dlambda
	// p_theta = dtheta/dlambda
	// p_phi = dphi/dlambda
};