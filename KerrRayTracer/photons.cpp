#include "photons.hpp"
#include "parameters.hpp"
#include "config.hpp"
#include "rendering_physics.hpp"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <omp.h>

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
	accumulatedColour.resize(i);
	transmittance.resize(i);
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
	constexpr int MAX_STEPS = 3;

	PhotonDerivative k2, k3, k4, k5, k6, k7;

	double next_r;
	double next_theta;
	for (int rejectedSteps = 0; rejectedSteps < MAX_STEPS; ++rejectedSteps) {
		//k1 = evaluate(r[i], theta[i], sign_r[i], sign_theta[i], xi[i], eta[i]);
		double r_stage = r[i] + (dlambda[i] / 5.0) * (dr[i] / 5.0);
		double theta_stage = theta[i] + (dlambda[i] / 5.0) * (dtheta[i] / 5.0);

		k2 = evaluate(r_stage, theta_stage, sign_r[i], sign_theta[i], xi[i], eta[i]);
		r_stage = r[i] + (dlambda[i] * 3.0 / 10.0) * (3.0 / 40.0 * dr[i] + 9.0 / 40.0 * k2.dr);
		theta_stage = theta[i] + (dlambda[i] * 3.0 / 10.0) * (3.0 / 40.0 * dtheta[i] + 9.0 / 40.0 * k2.dtheta);

		k3 = evaluate(r_stage, theta_stage, sign_r[i], sign_theta[i], xi[i], eta[i]);
		r_stage = r[i] + (dlambda[i] * 4.0 / 5.0) * (44.0 / 45.0 * dr[i] + -56.0 / 15.0 * k2.dr + 32.0 / 9.0 * k3.dr);
		theta_stage = theta[i] + (dlambda[i] * 4.0 / 5.0) * (44.0 / 45.0 * dtheta[i] + -56.0 / 15.0 * k2.dtheta + 32.0 / 9.0 * k3.dtheta);

		k4 = evaluate(r_stage, theta_stage, sign_r[i], sign_theta[i], xi[i], eta[i]);
		r_stage = r[i] + (dlambda[i] * 8.0 / 9.0) * (19372.0 / 6561.0 * dr[i] + -25360.0 / 2187.0 * k2.dr + 64448.0 / 6561.0 * k3.dr + -212.0 / 729.0 * k4.dr);
		theta_stage = theta[i] + (dlambda[i] * 8.0 / 9.0) * (19372.0 / 6561.0 * dtheta[i] + -25360.0 / 2187.0 * k2.dtheta + 64448.0 / 6561.0 * k3.dtheta + -212.0 / 729.0 * k4.dtheta);

		k5 = evaluate(r_stage, theta_stage, sign_r[i], sign_theta[i], xi[i], eta[i]);
		r_stage = r[i] + dlambda[i] * (9017.0 / 3168.0 * dr[i] + -355.0 / 33.0 * k2.dr + 46732.0 / 5247.0 * k3.dr + 49.0 / 176.0 * k4.dr + -5103.0 / 18656.0 * k5.dr);
		theta_stage = theta[i] + dlambda[i] * (9017.0 / 3168.0 * dtheta[i] + -355.0 / 33.0 * k2.dtheta + 46732.0 / 5247.0 * k3.dtheta + 49.0 / 176.0 * k4.dtheta + -5103.0 / 18656.0 * k5.dtheta);

		k6 = evaluate(r_stage, theta_stage, sign_r[i], sign_theta[i], xi[i], eta[i]);
	    next_r = r[i] + dlambda[i] * (35.0 / 384.0 * dr[i] + 500.0 / 1113.0 * k3.dr + 125.0 / 192.0 * k4.dr + -2187.0 / 6784.0 * k5.dr + 11.0 / 84.0 * k6.dr);
		next_theta = theta[i] + dlambda[i] * (35.0 / 384.0 * dtheta[i] + 500.0 / 1113.0 * k3.dtheta + 125.0 / 192.0 * k4.dtheta + -2187.0 / 6784.0 * k5.dtheta + 11.0 / 84.0 * k6.dtheta);

		k7 = evaluate(next_r, next_theta, sign_r[i], sign_theta[i], xi[i], eta[i]);
		double err_r = dlambda[i] * std::abs((32.0 / 384.0 - 5179.0 / 57600.0) * dr[i] + (500.0 / 1113.0 - 7371.0 / 16695.0) * k3.dr + (125.0 / 192.0 - 393.0 / 640.0) * k4.dr + (-2187.0 / 6784.0 + 92097.0 / 339200.0) * k5.dr + (11.0 / 84.0 - 187.0 / 2100.0) * k6.dr - k7.dr / 40.0);
		double err_theta = dlambda[i] * std::abs((32.0 / 384.0 - 5179.0 / 57600.0) * dtheta[i] + (500.0 / 1113.0 - 7371.0 / 16695.0) * k3.dtheta + (125.0 / 192.0 - 393.0 / 640.0) * k4.dtheta + (-2187.0 / 6784.0 + 92097.0 / 339200.0) * k5.dtheta + (11.0 / 84.0 - 187.0 / 2100.0) * k6.dtheta - k7.dtheta / 40.0);

		double norm_err = std::max(err_r / (atol + rtol * std::abs(r[i])), err_theta / (atol + rtol * std::abs(theta[i])));
		double scale = safety * std::pow(1.0 / std::max(norm_err, 1e-10), 0.2);
		scale = std::clamp(scale, min_scale, max_scale);

		dlambda[i] *= scale;
		if (norm_err <= 1.0 || rejectedSteps == MAX_STEPS - 1)
			break;
		//dlambda[i] *= scale;

		//if (norm_err <= 1.0 || rejectedSteps == MAX_STEPS - 1) {
		//	r[i] = std::max(next_r, 0.0);
		//	theta[i] = std::clamp(next_theta, 1e-14, PI - 1e-14);
		//	phi[i] += dlambda[i] * (35.0 / 384.0 * dphi[i] + 500.0 / 1113.0 * k3.dphi + 125.0 / 192.0 * k4.dphi + -2187.0 / 6784.0 * k5.dphi + 11.0 / 84.0 * k6.dphi);
		//	//phi[i] = std::clamp(phi[i], 0.0, TWO_PI); // do this to get the cool wavy glitched skyfield look
		//	dr[i] = k7.dr; 
		//	dtheta[i] = k7.dtheta;
		//	dphi[i] = k7.dphi;
		//	sign_r[i] *= (2 * (abs(dr[i]) >= 1e-14) - 1);
		//	sign_theta[i] *= (2 * (abs(dtheta[i]) >= 1e-14) - 1);
		//	r[i] += sign_r[i] * (abs(dr[i]) <= 1e-14) * dlambda[i];
		//	theta[i] += sign_theta[i] * (abs(dtheta[i]) <= 1e-14) * dlambda[i];
		//	break;
		//}
	}
	r[i] = std::max(next_r, 0.0);
	theta[i] = std::clamp(next_theta, 1e-14, PI - 1e-14);
	phi[i] += dlambda[i] * (35.0 / 384.0 * dphi[i] + 500.0 / 1113.0 * k3.dphi + 125.0 / 192.0 * k4.dphi + -2187.0 / 6784.0 * k5.dphi + 11.0 / 84.0 * k6.dphi);
	//phi[i] = std::clamp(phi[i], 0.0, TWO_PI); // do this to get the cool wavy glitched skyfield look
	dr[i] = k7.dr;
	dtheta[i] = k7.dtheta;
	dphi[i] = k7.dphi;
	sign_r[i] *= (2 * (abs(dr[i]) >= 1e-14) - 1);
	sign_theta[i] *= (2 * (abs(dtheta[i]) >= 1e-14) - 1);
	r[i] += sign_r[i] * (abs(dr[i]) <= 1e-14) * dlambda[i];
	theta[i] += sign_theta[i] * (abs(dtheta[i]) <= 1e-14) * dlambda[i];
}

void Photons::adaptiveRK4StepRay(size_t i) {
	PhotonDerivative k1, k2, k3, k4;
	//double dlambda = 0.1 * std::log10(r_horizon - 1.01);

	k1 = evaluate(r[i], theta[i], sign_r[i], sign_theta[i], xi[i], eta[i]);
	double delta = r[i] * r[i] - 2.0 * r[i] + a * a;
	double cosTh = cos(theta[i]);
	double sigma = r[i] * r[i] + a * a * cosTh * cosTh;
	double dlambda = std::min(0.04, 0.01 * std::max(delta, 0.01) / ((sigma + 1e-5) * (abs(k1.dr) + abs(k1.dtheta) + abs(k1.dphi) + 1.0)));

	double r_stage = r[i] + k1.dr * (dlambda / 2.0);
	double theta_stage = theta[i] + k1.dtheta * (dlambda / 2.0);

	k2 = evaluate(r_stage, theta_stage, sign_r[i], sign_theta[i], xi[i], eta[i]);
	r_stage = r[i] + k2.dr * (dlambda / 2.0);
	theta_stage = theta[i] + k2.dtheta * (dlambda / 2.0);

	k3 = evaluate(r_stage, theta_stage, sign_r[i], sign_theta[i], xi[i], eta[i]);
	r_stage = r[i] + k3.dr * dlambda;
	theta_stage = theta[i] + k3.dtheta * dlambda;

	k4 = evaluate(r_stage, theta_stage, sign_r[i], sign_theta[i], xi[i], eta[i]);

	double next_r = r[i] + dlambda / 6.0 * (k1.dr + 2.0 * k2.dr + 2.0 * k3.dr + k4.dr);
	double next_theta = theta[i] + dlambda / 6.0 * (k1.dtheta + 2.0 * k2.dtheta + 2.0 * k3.dtheta + k4.dtheta);

	r[i] = std::max(next_r, 0.0);
	theta[i] = std::clamp(next_theta, 1e-14, PI - 1e-14);
	phi[i] += dlambda / 6.0 * (k1.dphi + 2.0 * k2.dphi + 2.0 * k3.dphi + k4.dphi);

	sign_r[i] *= (2 * (abs(k4.dr) >= 1e-14) - 1);
	sign_theta[i] *= (2 * (abs(k4.dtheta) >= 1e-14) - 1);
	r[i] += sign_r[i] * (abs(k4.dr) <= 1e-14) * dlambda;
	theta[i] += sign_theta[i] * (abs(k4.dtheta) <= 1e-14) * dlambda;
}

void Photons::adaptiveRK5StepRay(size_t i) {
	PhotonDerivative k1, k2, k3, k4, k5, k6;

	k1 = evaluate(r[i], theta[i], sign_r[i], sign_theta[i], xi[i], eta[i]);
	double delta = r[i] * r[i] - 2.0 * r[i] + a * a;
	double cosTh = cos(theta[i]);
	double sigma = r[i] * r[i] + a * a * cosTh * cosTh;
	double dlambda = std::min(0.04, 0.01 * std::max(delta, 0.01) / ((sigma + 1e-5) * (abs(k1.dr) + abs(k1.dtheta) + abs(k1.dphi) + 1.0)));

	double r_stage = r[i] + (dlambda / 4.0) * (k1.dr /4.0);
	double theta_stage = theta[i] + (dlambda / 4.0) * (k1.dtheta / 4.0);

	k2 = evaluate(r_stage, theta_stage, sign_r[i], sign_theta[i], xi[i], eta[i]);
	r_stage = r[i] + (dlambda / 4.0) * (k1.dr / 8.0 + k2.dr / 8.0);
	theta_stage = theta[i] + (dlambda / 4.0) * (k1.dtheta / 8.0 + k2.dtheta / 8.0);

	k3 = evaluate(r_stage, theta_stage, sign_r[i], sign_theta[i], xi[i], eta[i]);
	r_stage = r[i] + (dlambda / 2.0) * (-k2.dr / 2.0);
	theta_stage = theta[i] + (dlambda / 2.0) * (-k2.dtheta / 2.0);

	k4 = evaluate(r_stage, theta_stage, sign_r[i], sign_theta[i], xi[i], eta[i]);
	r_stage = r[i] + (dlambda * 3.0 / 4.0) * (k1.dr * 3.0 / 16.0 + k4.dr * 9.0 / 16.0);
	theta_stage = theta[i] + (dlambda * 3.0 / 4.0) * (k1.dtheta * 3.0 / 16.0 + k4.dtheta * 9.0 / 16.0);

	k5 = evaluate(r_stage, theta_stage, sign_r[i], sign_theta[i], xi[i], eta[i]);
	r_stage = r[i] + dlambda * (k1.dr * -3.0 / 7.0 + k2.dr * 2.0 / 7.0 + k3.dr * 12.0 / 7.0 + k4.dr * -12.0 / 7.0 + k5.dr * 8.0 / 7.0);
	theta_stage = theta[i] + dlambda * (k1.dtheta * -3.0 / 7.0 + k2.dtheta * 2.0 / 7.0 + k3.dtheta * 12.0 / 7.0 + k4.dtheta * -12.0 / 7.0 + k5.dtheta * 8.0 / 7.0);

	k6 = evaluate(r_stage, theta_stage, sign_r[i], sign_theta[i], xi[i], eta[i]);
	double next_r = r[i] + (dlambda / 90.0) * (7.0 * k1.dr + 32.0 * k3.dr + 12.0 * k4.dr + 32.0 * k5.dr + 7.0 * k6.dr);
	double next_theta = theta[i] + (dlambda / 90.0) * (7.0 * k1.dtheta + 32.0 * k3.dtheta + 12.0 * k4.dtheta + 32.0 * k5.dtheta + 7.0 * k6.dtheta);

	r[i] = std::max(next_r, 0.0);
	theta[i] = std::clamp(next_theta, 1e-14, PI - 1e-14);
	phi[i] += (dlambda / 90.0) * (7.0 * k1.dphi + 32.0 * k3.dphi + 12.0 * k4.dphi + 32.0 * k5.dphi + 7.0 * k6.dphi);

	sign_r[i] *= (2 * (abs(k6.dr) >= 1e-14) - 1);
	sign_theta[i] *= (2 * (abs(k6.dtheta) >= 1e-14) - 1);
	r[i] += sign_r[i] * (abs(k6.dr) <= 1e-14) * dlambda;
	theta[i] += sign_theta[i] * (abs(k6.dtheta) <= 1e-14) * dlambda;
}

void Photons::traceAllRays() {
	constexpr int MAX_STEPS = 40000;

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
			if constexpr (ENABLE_OUTER_HORIZON) {
				if (r[i] <= r_horizon + 0.01) {
					state[i] = PhotonState::Captured;
					break;
				}
			}

			double oldTheta = theta[i];

			if constexpr (RKDP_INTEGRATION) {
				rkdpStepRay(i);
			}

			if constexpr (RK4_INTEGRATION) {
				adaptiveRK4StepRay(i);
			}

			if constexpr (RK5_INTEGRATION) {
				adaptiveRK5StepRay(i);
			}

			
			if constexpr (ENABLE_ACCRETION_DISK) {
				//if (intersectAccretionDisk(r[i], theta[i])) {
				if constexpr (ENABLE_OPAQUE_DISK) {
					if (intersectAccretionDisk(r[i], oldTheta, theta[i])) {
						state[i] = PhotonState::AccretionDiskHit;
						break;
						// accretion disk sampler algorithm
					}
				}
				if constexpr (!ENABLE_OPAQUE_DISK) {
					if (intersectAccretionDisk(r[i], oldTheta, theta[i])) {
						float localOpacity = std::clamp(diskDensity(r[i], phi[i], 0), 0.0f, 1.0f);
						float3 localEmission = rgbaToFloat3(temperatureToRGB(observedTemperature(r[i], xi[i])));
						localEmission = relativisticBeaming(r[i], xi[i], localEmission);
						(accumulatedColour[i])[0] += transmittance[i] * localEmission[0] * localOpacity;
						(accumulatedColour[i])[1] += transmittance[i] * localEmission[1] * localOpacity;
						(accumulatedColour[i])[2] += transmittance[i] * localEmission[2] * localOpacity;
						transmittance[i] *= (1.0f - localOpacity);
					}
				}
			}
		}
	}
}