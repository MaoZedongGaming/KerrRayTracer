#include "photons.hpp"
#include "parameters.hpp"
#include "config.hpp"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <omp.h>

double accretionThickness(double r) {
	constexpr double H0 = 0.02;
	return H0 * pow((r / r_ISCO), 1.2);
}

constexpr bool crossedEquatorialPlane(double theta_old, double theta_new) {
	//return (std::cos(theta_old) * std::cos(theta_new) <= 0.0);
	return (std::min(theta_old, theta_new) <= PI_2 && PI_2 <= std::max(theta_old, theta_new));
}

//bool intersectAccretionDisk(double r, double theta) {
//	return ((abs(r * cos(theta)) <= accretionThickness(r) + 1e-3) && (r_ISCO <= r && r <= r_acceretion));
//}

bool intersectAccretionDisk(double r, double theta_old, double theta_new) {
	return crossedEquatorialPlane(theta_old, theta_new) && (r_ISCO < r && r <= r_acceretion);
}

void Photons::resize(size_t i) {
	r.resize(i);
	theta.resize(i);
	phi.resize(i);
	dr.resize(i);
	dtheta.resize(i);
	dphi.resize(i);
	xi.resize(i);
	eta.resize(i);
	sign_r.resize(i);
	sign_theta.resize(i);
	state.resize(i);
	dlambda.resize(i);
	count = i;
}

PhotonDerivative evaluate(double r, double theta, float sign_r, float sign_theta, double xi, double eta) {
	double r2 = r * r;
	double a2 = a * a;
	double cosTh = cos(theta);
	double cos2 = cosTh * cosTh;
	double sinTh = sin(theta);
	double sin2 = sinTh * sinTh;
	double delta = r2 - 2.0 * r + a2;
	double P = r2 + a2 - a * xi;
	double R = P * P - delta * (eta + (xi - a) * (xi - a));
	double Theta = eta + a2 * cos2 - xi * xi * cos2 / sin2;

	double dr = sign_r * sqrt(std::max(R, 0.0));
	double dtheta = sign_theta * sqrt(std::max(Theta, 0.0));
	double dphi = -(a - xi / sin2) + P * a / delta;

	return PhotonDerivative{ dr, dtheta, dphi };
}

void Photons::rkdpStepRay(size_t i) {
	constexpr double rtol = 1e-6;
	constexpr double atol = 1e-6;
	constexpr double safety = 0.9;
	constexpr double min_scale = 0.2;
	constexpr double max_scale = 10.0;
	constexpr int MAX_STEPS = 5;

	PhotonDerivative k2, k3, k4, k5, k6, k7;
	for (int rejectedSteps = 0; rejectedSteps < MAX_STEPS; ++rejectedSteps) {
		//k1 = evaluate(r[i], theta[i], sign_r[i], sign_theta[i], xi[i], eta[i]);
		double r_stage = r[i] + dlambda[i] * (dr[i] / 5.0);
		double theta_stage = theta[i] + dlambda[i] * (dtheta[i] / 5.0);

		k2 = evaluate(r_stage, theta_stage, sign_r[i], sign_theta[i], xi[i], eta[i]);
		r_stage = r[i] + dlambda[i] * (3.0 / 40.0 * dr[i] + 9.0 / 40.0 * k2.dr);
		theta_stage = theta[i] + dlambda[i] * (3.0 / 40.0 * dtheta[i] + 9.0 / 40.0 * k2.dtheta);

		k3 = evaluate(r_stage, theta_stage, sign_r[i], sign_theta[i], xi[i], eta[i]);
		r_stage = r[i] + dlambda[i] * (44.0 / 45.0 * dr[i] + -56.0 / 15.0 * k2.dr + 32.0 / 9.0 * k3.dr);
		theta_stage = theta[i] + dlambda[i] * (44.0 / 45.0 * dtheta[i] + -56.0 / 15.0 * k2.dtheta + 32.0 / 9.0 * k3.dtheta);

		k4 = evaluate(r_stage, theta_stage, sign_r[i], sign_theta[i], xi[i], eta[i]);
		r_stage = r[i] + dlambda[i] * (19372.0 / 6561.0 * dr[i] + -25360.0 / 2187.0 * k2.dr + 64448.0 / 6561.0 * k3.dr + -212.0 / 729.0 * k4.dr);
		theta_stage = theta[i] + dlambda[i] * (19372.0 / 6561.0 * dtheta[i] + -25360.0 / 2187.0 * k2.dtheta + 64448.0 / 6561.0 * k3.dtheta + -212.0 / 729.0 * k4.dtheta);

		k5 = evaluate(r_stage, theta_stage, sign_r[i], sign_theta[i], xi[i], eta[i]);
		r_stage = r[i] + dlambda[i] * (9017.0 / 3168.0 * dr[i] + -355.0 / 33.0 * k2.dr + 46732.0 / 5247.0 * k3.dr + 49.0 / 176.0 * k4.dr + -5103.0 / 18656.0 * k5.dr);
		theta_stage = theta[i] + dlambda[i] * (9017.0 / 3168.0 * dtheta[i] + -355.0 / 33.0 * k2.dtheta + 46732.0 / 5247.0 * k3.dtheta + 49.0 / 176.0 * k4.dtheta + -5103.0 / 18656.0 * k5.dtheta);

		k6 = evaluate(r_stage, theta_stage, sign_r[i], sign_theta[i], xi[i], eta[i]);
		double next_r = r[i] + dlambda[i] * (35.0 / 384.0 * dr[i] + 500.0 / 1113.0 * k3.dr + 125.0 / 192.0 * k4.dr + -2187.0 / 6784.0 * k5.dr + 11.0 / 84.0 * k6.dr);
		double next_theta = theta[i] + dlambda[i] * (35.0 / 384.0 * dtheta[i] + 500.0 / 1113.0 * k3.dtheta + 125.0 / 192.0 * k4.dtheta + -2187.0 / 6784.0 * k5.dtheta + 11.0 / 84.0 * k6.dtheta);

		k7 = evaluate(next_r, next_theta, sign_r[i], sign_theta[i], xi[i], eta[i]);
		double err_r = dlambda[i] * std::abs((32.0 / 384.0 - 5179.0 / 57600.0) * dr[i] + (500.0 / 1113.0 - 7371.0 / 16695.0) * k3.dr + (125.0 / 192.0 - 393.0 / 640.0) * k4.dr + (-2187.0 / 6784.0 + 92097.0 / 339200.0) * k5.dr + (11.0 / 84.0 - 187.0 / 2100.0) * k6.dr - k7.dr / 40.0);
		double err_theta = dlambda[i] * std::abs((32.0 / 384.0 - 5179.0 / 57600.0) * dtheta[i] + (500.0 / 1113.0 - 7371.0 / 16695.0) * k3.dtheta + (125.0 / 192.0 - 393.0 / 640.0) * k4.dtheta + (-2187.0 / 6784.0 + 92097.0 / 339200.0) * k5.dtheta + (11.0 / 84.0 - 187.0 / 2100.0) * k6.dtheta - k7.dtheta / 40.0);

		double norm_err = std::max(err_r / (atol + rtol * std::abs(r[i])), err_theta / (atol + rtol * std::abs(theta[i])));
		double scale = safety * std::pow(1.0 / std::max(norm_err, 1e-10), 0.2);
		scale = std::clamp(scale, min_scale, max_scale);

		dlambda[i] *= scale;

		if (norm_err <= 1.0 || rejectedSteps == MAX_STEPS - 1) {
			r[i] = std::max(next_r, 0.0);
			theta[i] = std::clamp(next_theta, 1e-14, PI - 1e-14);
			phi[i] += dlambda[i] * (35.0 / 384.0 * dphi[i] + 500.0 / 1113.0 * k3.dphi + 125.0 / 192.0 * k4.dphi + -2187.0 / 6784.0 * k5.dphi + 11.0 / 84.0 * k6.dphi);
			//phi[i] = std::clamp(phi[i], 0.0, TWO_PI); // do this to get the cool wavy glitched skyfield look
			dr[i] = k7.dr; 
			dtheta[i] = k7.dtheta;
			dphi[i] = k7.dphi;
			dr[i] *= (2 * (abs(dr[i]) >= 1e-14) - 1);
			dtheta[i] *= (2 * (abs(dtheta[i]) >= 1e-14) - 1);
			sign_r[i] *= (2 * (abs(dr[i]) >= 1e-14) - 1);
			sign_theta[i] *= (2 * (abs(dtheta[i]) >= 1e-14) - 1);
			r[i] += sign_r[i] * (abs(dr[i]) <= 1e-14) * dlambda[i];
			theta[i] += sign_theta[i] * (abs(dtheta[i]) <= 1e-14) * dlambda[i];
			break;
		}
	}
}

void Photons::traceAllRays() {
	constexpr int MAX_STEPS = 30000;

	#pragma omp parallel for schedule(dynamic, 16)
	for (int i = 0; i < count; ++i) {
		for (size_t steps = 0; steps < MAX_STEPS; ++steps) {
			//constexpr size_t debugIdx = 599 + (600 * 200);
			//if (i == debugIdx) {
			//	std::cout << "r theta phi = (" << r[debugIdx] << ", " << theta[debugIdx] << ", " << phi[debugIdx] << ") \n";
			//	//std::cout << "phi = " << phi[304] << " dr = " << dphi[304] << "\n";
			//	std::cout << "dr dtheta dphi = (" << dr[debugIdx] << ", " << dtheta[debugIdx] << ", " << dphi[debugIdx] << ") \n";
			//}
			if (r[i] >= r_sky) {
				state[i] = PhotonState::Escaped;
				break;
			}
			if (r[i] <= r_horizon + 0.01) {
				state[i] = PhotonState::Captured;
				break;
			}

			double oldTheta = theta[i];
			rkdpStepRay(i);
			
			if constexpr (ENABLE_ACCRETION_DISK) {
				//if (intersectAccretionDisk(r[i], theta[i])) {
				if (intersectAccretionDisk(r[i], oldTheta, theta[i])) {
					state[i] = PhotonState::AccretionDiskHit;
					break;
					// accretion disk sampler algorithm
				}
			}
		}
	}
}