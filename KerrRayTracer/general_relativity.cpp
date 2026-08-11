#include "general_relativity.hpp"

template <typename T>
T metricAt(size_t mu, size_t nu, Vector4d const& pos, T mass, T spin) {
	double r = pos[1];
	double theta = pos[2];
	double sigma_val = sigma(r, mass, spin);
	double delta_val = delta(r, mass, spin);
	if (mu == 0 && nu == 0) {
		return -(1 - 2 * mass * r / sigma_val);
	}
	else if ((mu == 0 && nu == 3) || (mu == 3 && nu == 0)) {
		return -2 * mass * r * spin * std::sin(theta) * std::sin(theta) / sigma_val;
	}
	else if (mu == 1 && nu == 1) {
		return sigma_val / delta_val;
	}
	else if (mu == 2 && nu == 2) {
		return sigma_val;
	}
	else if (mu == 3 && nu == 3) {
		return (r * r + spin * spin + 2 * mass * r * spin * spin * std::sin(theta) * std::sin(theta) / sigma_val) * std::sin(theta) * std::sin(theta);
	}
	else {
		return T(0);
	}
}