#include "general_relativity.hpp"
#include "length_scales.hpp"

double metricAt(size_t mu, size_t nu, double r, double theta) {
	double sigma_val = sigma(r, theta);
	double delta_val = delta(r);
	if (mu == 0 && nu == 0) {
		return -(1 - 2 * r / sigma_val);
	}
	else if ((mu == 0 && nu == 3) || (mu == 3 && nu == 0)) {
		return -4 * r * a * std::sin(theta) * std::sin(theta) / sigma_val;
	}
	else if (mu == 1 && nu == 1) {
		return sigma_val / delta_val;
	}
	else if (mu == 2 && nu == 2) {
		return sigma_val;
	}
	else if (mu == 3 && nu == 3) {
		return (r * r + a * a + 2 * r * a * a * std::sin(theta) * std::sin(theta) / sigma_val) * std::sin(theta) * std::sin(theta);
	}
	else {
		return 0.0;
	}
}

//double metricDerivativeAt(size_t alpha, size_t mu, size_t nu, Vector4d const& pos, double mass, double spin) {
//	if (mu != nu && (mu != 0 && nu != 3) && (mu != 3 && nu != 0))
//		return 0.0;
//	Vector4X<Dual<double>> dualPos;
//	for (size_t i = 0; i < 4; ++i)
//		dualPos[i] = Dual<double>(pos[i], (i == alpha) ? 1.0 : 0.0);
//	return metricAt(mu, nu, dualPos, mass, spin).derivative;
//}

double inverseMetricAt(size_t mu, size_t nu, double r, double theta) {
	double sigma_val = sigma(r, theta);
	double delta_val = delta(r);
	if (mu == 0 && nu == 0) {
		return -(r * r + a * a + 2 * r * a * a * std::sin(theta) * std::sin(theta) / sigma_val) / (sigma_val * delta_val);
	}
	else if ((mu == 0 && nu == 3) || (mu == 3 && nu == 0)) {
		return -4 * r * a / (sigma_val * delta_val);
	}
	else if (mu == 1 && nu == 1) {
		return delta_val / sigma_val;
	}
	else if (mu == 2 && nu == 2) {
		return 1.0 / sigma_val;
	}
	else if (mu == 3 && nu == 3) {
		return (1 - 2 * r / sigma_val) / (sigma_val * std::sin(theta) * std::sin(theta));
	}
	else {
		return 0.0;
	}
}

//double christoffelSymbol(size_t lambda, size_t mu, size_t nu, Vector4d const& pos, double mass, double spin) {
//	double res = 0.0;
//	for (size_t alpha = 0; alpha < 4; ++alpha) {  //sum over alpha because non diagonal elements of the metric are non-zero
//		double g_inv = inverseMetricAt(lambda, alpha, pos, mass, spin);
//		double dg_mu_nu = metricDerivativeAt(alpha, mu, nu, pos, mass, spin);
//		double dg_mu_alpha = metricDerivativeAt(mu, alpha, nu, pos, mass, spin);
//		double dg_nu_alpha = metricDerivativeAt(nu, alpha, mu, pos, mass, spin);
//		res += g_inv * (dg_mu_nu + dg_nu_alpha - dg_mu_alpha);
//	}
//	return 0.5 * res;
//}