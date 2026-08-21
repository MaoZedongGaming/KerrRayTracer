#include "photons.hpp"
#include "parameters.hpp"
#include "length_scales.hpp"
#include <cmath>
#include <iostream>
#include <algorithm>
#include "omp.h"


constexpr double PI_2 = 1.57079632679;

template <typename T>
constexpr int sgn(T x) {
	return ((x > T(0)) - (x < T(0)));
}

template <>
constexpr int sgn(double x) {
	return ((x > 1e-14) - (1e-14 > x));
}

double accretionThickness(double r) {
	constexpr double H0 = 0.05;
	return H0 * pow((r / r_ISCO), 1.2);
}

bool intersectAccretionDisk(double r, double theta) {
	return ((abs(r * cos(theta)) <= accretionThickness(r)) && (r_ISCO <= r || r <= r_acceretion));
}

void Photons::reserve(size_t i) {
	r.reserve(i);
	theta.reserve(i);
	phi.reserve(i);
	xi.reserve(i);
	eta.reserve(i);
	sign_r.reserve(i);
	sign_theta.reserve(i);
	state.reserve(i);
	//activeIndices.reserve(i);
	dlambda.reserve(i);
}

void Photons::resize(size_t i) {
	r.resize(i);
	theta.resize(i);
	phi.resize(i);
	xi.resize(i);
	eta.resize(i);
	sign_r.resize(i);
	sign_theta.resize(i);
	state.resize(i);
	//activeIndices.resize(i);
	dlambda.resize(i);
}

void Photons::clear() {
	r.clear();
	theta.clear();
	phi.clear();
	xi.clear();
	eta.clear();
	sign_r.clear();
	sign_theta.clear();
	state.clear();
	//activeIndices.clear();
	dlambda.clear();
	count = 0;
}

struct PhotonDerivative {
	//double dt;
	double dr;
	double dtheta;
	double dphi;
};

PhotonDerivative evaluate(double r, double theta, float sign_r, float sign_theta, double xi, double eta) {
	double P = r * r + a * a - a * xi;
	double R = P * P - delta(r) * (eta + (xi - a) * (xi - a));
	double Theta = eta + a * a * cos(theta) * cos(theta) - xi * xi / (tan(theta) * tan(theta));

	//double dt = -a * (a * sin(theta) * sin(theta) - xi) + (r * r + a * a) * P / delta(r);
	double dr = sign_r * sqrt(std::max(0.0, R));
	double dtheta = sign_theta * sqrt(std::max(0.0, Theta));
	double dphi = -(a - xi / (sin(theta) * sin(theta))) + P * a / delta(r);

	return PhotonDerivative{ /*dt,*/ dr, dtheta, dphi };
}

void Photons::rkdpStepRay(size_t i) {
	constexpr double rtol = 1e-6;
	constexpr double atol = 1e-6;
	constexpr double safety = 0.9;
	constexpr double min_scale = 0.2;
	constexpr double max_scale = 10.0;
	constexpr int MAX_STEPS = 50;

	PhotonDerivative k1, k2, k3, k4, k5, k6, k7;
	for (int rejectedSteps = 0; rejectedSteps < MAX_STEPS; ++rejectedSteps) {
		k1 = evaluate(r[i], theta[i], sign_r[i], sign_theta[i], xi[i], eta[i]);
		double r_stage = r[i] + dlambda[i] * (k1.dr / 5.0);
		double theta_stage = theta[i] + dlambda[i] * (k1.dtheta / 5.0);

		k2 = evaluate(r_stage, theta_stage, sign_r[i], sign_theta[i], xi[i], eta[i]);
		r_stage = r[i] + dlambda[i] * (3.0 / 40.0 * k1.dr + 9.0 / 40.0 * k2.dr);
		theta_stage = theta[i] + dlambda[i] * (3.0 / 40.0 * k1.dtheta + 9.0 / 40.0 * k2.dtheta);

		k3 = evaluate(r_stage, theta_stage, sign_r[i], sign_theta[i], xi[i], eta[i]);
		r_stage = r[i] + dlambda[i] * (44.0 / 45.0 * k1.dr + -56.0 / 15.0 * k2.dr + 32.0 / 9.0 * k3.dr);
		theta_stage = theta[i] + dlambda[i] * (44.0 / 45.0 * k1.dtheta + -56.0 / 15.0 * k2.dtheta + 32.0 / 9.0 * k3.dtheta);

		k4 = evaluate(r_stage, theta_stage, sign_r[i], sign_theta[i], xi[i], eta[i]);
		r_stage = r[i] + dlambda[i] * (19372.0 / 6561.0 * k1.dr + -25360.0 / 2187.0 * k2.dr + 64448.0 / 6561.0 * k3.dr + -212.0 / 729.0 * k4.dr);
		theta_stage = theta[i] + dlambda[i] * (19372.0 / 6561.0 * k1.dtheta + -25360.0 / 2187.0 * k2.dtheta + 64448.0 / 6561.0 * k3.dtheta + -212.0 / 729.0 * k4.dtheta);

		k5 = evaluate(r_stage, theta_stage, sign_r[i], sign_theta[i], xi[i], eta[i]);
		r_stage = r[i] + dlambda[i] * (9017.0 / 3168.0 * k1.dr + -355.0 / 33.0 * k2.dr + 46732.0 / 5247.0 * k3.dr + 49.0 / 176.0 * k4.dr + -5103.0 / 18656.0 * k5.dr);
		theta_stage = theta[i] + dlambda[i] * (9017.0 / 3168.0 * k1.dtheta + -355.0 / 33.0 * k2.dtheta + 46732.0 / 5247.0 * k3.dtheta + 49.0 / 176.0 * k4.dtheta + -5103.0 / 18656.0 * k5.dtheta);

		k6 = evaluate(r_stage, theta_stage, sign_r[i], sign_theta[i], xi[i], eta[i]);
		double next_r = r[i] + dlambda[i] * (35.0 / 384.0 * k1.dr + 500.0 / 1113.0 * k3.dr + 125.0 / 192.0 * k4.dr + -2187.0 / 6784.0 * k5.dr + 11.0 / 84.0 * k6.dr);
		double next_theta = theta[i] + dlambda[i] * (35.0 / 384.0 * k1.dtheta + 500.0 / 1113.0 * k3.dtheta + 125.0 / 192.0 * k4.dtheta + -2187.0 / 6784.0 * k5.dtheta + 11.0 / 84.0 * k6.dtheta);

		k7 = evaluate(next_r, next_theta, sign_r[i], sign_theta[i], xi[i], eta[i]);
		double err_r = dlambda[i] * std::abs((32.0 / 384.0 - 5179.0 / 57600.0) * k1.dr + (500.0 / 1113.0 - 7371.0 / 16695.0) * k3.dr + (125.0 / 192.0 - 393.0 / 640.0) * k4.dr + (-2187.0 / 6784.0 + 92097.0 / 339200.0) * k5.dr + (11.0 / 84.0 - 187.0 / 2100.0) * k6.dr - k7.dr / 40.0);
		double err_theta = dlambda[i] * std::abs((32.0 / 384.0 - 5179.0 / 57600.0) * k1.dtheta + (500.0 / 1113.0 - 7371.0 / 16695.0) * k3.dtheta + (125.0 / 192.0 - 393.0 / 640.0) * k4.dtheta + (-2187.0 / 6784.0 + 92097.0 / 339200.0) * k5.dtheta + (11.0 / 84.0 - 187.0 / 2100.0) * k6.dtheta - k7.dtheta / 40.0);

		double norm_err = std::max(err_r / (atol + rtol * std::abs(r[i])), err_theta / (atol + rtol * std::abs(theta[i])));
		double scale = safety * std::pow(1.0 / std::max(norm_err, 1e-10), 0.2);
		scale = std::clamp(scale, min_scale, max_scale);

		dlambda[i] *= scale;

		if (norm_err <= 1.0 || rejectedSteps == MAX_STEPS) {
			//t[i] += dlambda[i] * (35.0 / 384.0 * k1.dt + 500.0 / 1113.0 * k3.dt + 125.0 / 192.0 * k4.dt + -2187.0 / 6784.0 * k5.dt + 11.0 / 84.0 * k6.dt);
			r[i] = next_r;
			theta[i] = next_theta;
			phi[i] += dlambda[i] * (35.0 / 384.0 * k1.dphi + 500.0 / 1113.0 * k3.dphi + 125.0 / 192.0 * k4.dphi + -2187.0 / 6784.0 * k5.dphi + 11.0 / 84.0 * k6.dphi);
			double P_next = r[i] * r[i] + a * a - a * xi[i];
			double R_next = P_next * P_next - delta(r[i]) * (eta[i] + (xi[i] - a) * (xi[i] - a));
			double Theta_next = eta[i] + a * a * cos(theta[i]) * cos(theta[i]) - xi[i] * xi[i] / (tan(theta[i]) * tan(theta[i]));
			sign_r[i] *= sgn(R_next);
			sign_theta[i] *= sgn(Theta_next);
			break;
		}
	}
}

void Photons::traceAllRays() {
	constexpr int MAX_STEPS = 1000;
	
	#pragma omp parallel for schedule(dynamic, 16)
	for (int i = 0; i < count; ++i) {
		for (size_t steps = 0; steps < MAX_STEPS; ++steps) {

			if (state[i] != PhotonState::Active)
				break;

			rkdpStepRay(i);

			if (r[i] >= r_sky) {
				state[i] = PhotonState::Escaped;
				//std::cout << "escaped! \n";
				// sample from skybox texture
			}
			if (r[i] <= r_horizon + 1e-4) {
				state[i] = PhotonState::Captured;
				//std::cout << "captured! \n";
				// black pixel
			}
			if (intersectAccretionDisk(r[i], theta[i])) {
				state[i] = PhotonState::AccretionDiskHit;
				// accretion disk sampler algorithm
			}
		}
	}
}