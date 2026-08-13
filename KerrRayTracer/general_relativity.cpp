#include "general_relativity.hpp"
#include "length_scales.hpp"

constexpr double g_tt(double r, double theta) {
	return -(1 - 2 * r / sigma(r, theta));
}

constexpr double g_rr(double r, double theta) {
	return sigma(r, theta)/delta(r);
}

constexpr double g_thth(double r, double theta) {
	return sigma(r, theta);
}

constexpr double g_phiphi(double r, double theta) {
	return (r * r + a * a + 2 * r * a * a * std::sin(theta) * std::sin(theta) / sigma(r, theta)) * std::sin(theta) * std::sin(theta);
}

constexpr double g_tphi(double r, double theta) {
	return -4 * r * a * std::sin(theta) * std::sin(theta) / sigma(r, theta);
}

