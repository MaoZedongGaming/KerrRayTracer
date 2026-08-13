#include "relativistic_camera.hpp"
#include "general_relativity.hpp"
#include "length_scales.hpp"
#include "photons.hpp"
#include "maths.hpp"
#include <cmath>

Vector4d lowerContravariant(Vector4d const& v, double r, double theta) {
	Vector4d res;
	res[0] = g_tt(r, theta) * v[0] + g_tphi(r, theta) * v[3];
	res[1] = g_rr(r, theta) * v[1];
	res[2] = g_thth(r, theta) * v[2];
	res[3] = g_tphi(r, theta) * v[0] + g_phiphi(r, theta) * v[3];
	return res;
}

double length2(Vector4d const& v, double r, double theta) {
	Vector4d v_lower = lowerContravariant(v, r, theta);
	double res = 0.0;
	for (size_t i = 0; i < 4; ++i)
		res += v_lower[i] * v[i];
	return res;
}

RelativisticCamera::RelativisticCamera(int w, int h) : width(w), height(h) {
	photons.reserve(w * h);
}

void RelativisticCamera::initTetrad() {
	//ZAMO frame 
	double r = position[1];
	double theta = position[2];

	
	double omega = -g_tphi(r, theta) / g_phiphi(r, theta);
	double alpha = std::sqrt((g_tphi(r, theta) * g_tphi(r, theta) / g_phiphi(r, theta)) - g_tt(r, theta));

	frame.e0 = Vector4d({ 1.0 / alpha, 0.0, 0.0, omega / alpha });
	Vector4d e1 = Vector4d({ 0.0, 1.0 / std::sqrt(g_rr(r, theta)), 0.0, 0.0 });
	Vector4d e2 = Vector4d({ 0.0, 0.0, 1.0 / std::sqrt(g_thth(r, theta)), 0.0 });
	Vector4d e3 = Vector4d({ 0.0, 0.0, 0.0, 1.0 / std::sqrt(g_phiphi(r, theta)) });
	double proj = g_tphi(r, theta) * e3[3] * frame.e0[0] + g_phiphi(r, theta) * e3[3] * frame.e0[3];
	e3 = e3 + frame.e0 * proj;
	e3 /= length2(e3, r, theta);

	frame.e1 = (cos(yaw) * cos(pitch)) * e1 + (cos(yaw) * sin(pitch) * sin(roll) - sin(yaw) * cos(roll)) * e2 + (cos(yaw) * sin(pitch) * cos(roll) + sin(yaw) * sin(roll)) * e3;
	frame.e2 = (sin(yaw) * cos(pitch)) * e1 + (sin(yaw) * sin(pitch) * sin(roll) - cos(yaw) * cos(roll)) * e2 + (sin(yaw) * sin(pitch) * cos(roll) + cos(yaw) * sin(roll)) * e3;
	frame.e3 = -sin(pitch) * e1 + (cos(pitch) * sin(roll)) * e2 + (cos(pitch) * cos(roll)) * e3;
}


void RelativisticCamera::generatePhotons() {
	double r = position[1];
	double theta = position[2];
	for (size_t j = 0; j < height; ++j) {
		for (size_t i = 0; i < width; ++i) {
			photons.t.push_back(position[0]);
			photons.r.push_back(position[1]);
			photons.theta.push_back(position[2]);
			photons.phi.push_back(position[3]);
			photons.activeIndices.push_back(i + width * j);
			photons.state.push_back(PhotonState::Active);
			photons.sign_r.push_back(1.0f);
			photons.sign_theta.push_back(1.0f);

			double screenX = (2.0 * i / (double)width - 1.0) * tan(fov / 2.0) * (width / (double)height);
			double screenY = (1.0 - 2.0 * j / (double)height) * tan(fov / 2.0);

			// momenta in the camera's tetrad frame, p_t = E = 1.0
			double p_1 = 1.0 / sqrt(1.0 + screenX * screenX + screenY * screenY); 
			double p_2 = screenY * p_1;
			double p_3 = screenX * p_1;

			// momentum projected onto global coordinates
			Vector4d p = frame.e0 + p_1 * frame.e1 + p_2 * frame.e2 + p_3 * frame.e3;

			// conserved covariant constants
			double E = p[0] * g_tt(r, theta) + p[3] * g_tphi(r, theta);
			double L_z = p[3] * g_phiphi(r, theta) + p[0] * g_tphi(r, theta);
			double Q = p[2] * p[2] * g_thth(r, theta) * g_thth(r, theta) + cos(theta) * cos(theta) * (a * a * (- E * E) + (L_z * L_z) / (sin(theta) * sin(theta)));

			photons.xi.push_back(L_z / E);
			photons.eta.push_back(Q / (E * E));
		}
	}
}