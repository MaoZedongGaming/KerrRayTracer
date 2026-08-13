#include "photons.hpp"

void Photons::reserve(size_t i) {
	t.reserve(i);
	r.reserve(i);
	theta.reserve(i);
	phi.reserve(i);
	xi.reserve(i);
	eta.reserve(i);
	sign_r.reserve(i);
	sign_theta.reserve(i);
	state.reserve(i);
	activeIndices.reserve(i);
}