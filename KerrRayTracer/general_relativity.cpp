#include "general_relativity.hpp"
#include "length_scales.hpp"
#include "parameters.hpp"

double g_tt(double r, double theta) {
	return -(1 - 2 * r / sigma(r, theta));
}

double g_rr(double r, double theta) {
	return sigma(r, theta)/delta(r);
}

double g_thth(double r, double theta) {
	return sigma(r, theta);
}

double g_phiphi(double r, double theta) {
	return (r * r + a * a + 2 * r * a * a * std::sin(theta) * std::sin(theta) / sigma(r, theta)) * std::sin(theta) * std::sin(theta);
}

double g_tphi(double r, double theta) {
	return -4 * r * a * std::sin(theta) * std::sin(theta) / sigma(r, theta);
}

