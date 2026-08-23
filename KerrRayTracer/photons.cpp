#include "photons.hpp"
#include "parameters.hpp"
#include "length_scales.hpp"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <omp.h>

double accretionThickness(double r) {
	constexpr double H0 = 1.0;
	return H0 * pow((r / r_ISCO), 1.2);
}

bool crossedEquatorialPlane(double theta_old, double theta_new) {
	//return (std::cos(theta_old) * std::cos(theta_new) <= 0.0);
	return ((theta_old <= PI_2 && PI_2 <= theta_new) || ((theta_new <= PI_2 && PI_2 <= theta_new)));
}

bool intersectAccretionDisk(double r, double theta) {
	return ((abs(r * cos(theta)) <= accretionThickness(r) + 1e-3) && (r_ISCO <= r && r <= r_acceretion));
}

//void Photons::reserve(size_t i) {
//	r.reserve(i);
//	theta.reserve(i);
//	phi.reserve(i);
//	xi.reserve(i);
//	eta.reserve(i);
//	sign_r.reserve(i);
//	sign_theta.reserve(i);
//	state.reserve(i);
//	//activeIndices.reserve(i);
//	dlambda.reserve(i);
//}

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

//void Photons::clear() {
//	r.clear();
//	theta.clear();
//	phi.clear();
//	xi.clear();
//	eta.clear();
//	sign_r.clear();
//	sign_theta.clear();
//	state.clear();
//	//activeIndices.clear();
//	dlambda.clear();
//	count = 0;
//}

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
			//std::cout << "finished RKDP step in " << rejectedSteps << " steps, with norm_err = " << norm_err << " \n";
			r[i] = next_r;
			bool thetaOOB = next_theta <= 1e-14 || next_theta >= PI - 1e-14;
			phi[i] += thetaOOB * PI;
			theta[i] = (2 * !thetaOOB - 1) * next_theta + (next_theta >= PI - 1e-14) * TWO_PI;
			phi[i] += dlambda[i] * (35.0 / 384.0 * dphi[i] + 500.0 / 1113.0 * k3.dphi + 125.0 / 192.0 * k4.dphi + -2187.0 / 6784.0 * k5.dphi + 11.0 / 84.0 * k6.dphi);
			phi[i] = std::fmod(phi[i], TWO_PI);
			phi[i] += std::signbit(phi[i]) * TWO_PI;
			dr[i] = k7.dr; 
			dtheta[i] = k7.dtheta;
			dphi[i] = k7.dphi;
			sign_r[i] *= (2 * (abs(dr[i]) > 1e-14) - 1);
			sign_theta[i] *= (2 * (abs(dtheta[i]) > 1e-14) - 1);
			break;
		}
	}
}

void Photons::traceAllRays() {
	constexpr int MAX_STEPS = 20000;

	#pragma omp parallel for schedule(dynamic, 16)
	for (int i = 0; i < count; ++i) {
		for (size_t steps = 0; steps < MAX_STEPS; ++steps) {
			if (r[i] >= r_sky) {
				state[i] = PhotonState::Escaped;
				break;
				//std::cout << "escaped! \n";
				// sample from skybox texture, this works
			}
			if (r[i] <= r_horizon + 0.05) {
				state[i] = PhotonState::Captured;
				break;
				//std::cout << "captured! \n";
				// black pixel
			}

			double oldTheta = theta[i];
			rkdpStepRay(i);
			
			//if (intersectAccretionDisk(r[i], theta[i])) {
			if (crossedEquatorialPlane(oldTheta, theta[i]) && (r_ISCO <= r[i] && r[i] <= r_acceretion)) {
				state[i] = PhotonState::AccretionDiskHit;
				break;
				//std::cout << "Hit AccretionDisk! \n";
				// accretion disk sampler algorithm
			}
		}
	}
}