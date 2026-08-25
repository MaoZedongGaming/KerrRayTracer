#include "relativistic_camera.hpp"
#include "general_relativity.hpp"
#include "parameters.hpp"
#include "photons.hpp"
#include "maths.hpp"
#include <algorithm>
#include <iostream>
#include <omp.h>
#include <cmath>

RelativisticCamera::RelativisticCamera(size_t w, size_t h) : width(w), height(h) {
	photons.resize(w * h);
	pixelBuffer.resize(w * h);
}

void RelativisticCamera::setPosition(double r, double theta) {
	position[0] = 0;
	position[1] = r;
	position[2] = std::clamp(theta, 0.01, PI - 0.01);
	position[3] = 0;

	velocity[0] = 1.0;
	velocity[1] = 0;
	velocity[2] = 0;
	velocity[3] = 0;
}

void RelativisticCamera::initTetrad() {
	//ZAMO frame 
	double r = position[1];
	double theta = position[2];

	double g_tp = g_tphi(r, theta);
	double g_pp = g_phiphi(r, theta);
	double g_t = g_tt(r, theta);
	
	double omega = -g_tp / g_pp;
	double alpha = std::sqrt((g_tp * g_tp / g_pp) - g_t);

	frame.e0 = Vector4d({ 1.0 / alpha, 0.0, 0.0, omega / alpha });
	Vector4d e1 = Vector4d({ 0.0, 1.0 / std::sqrt(g_rr(r, theta)), 0.0, 0.0 });
	Vector4d e2 = Vector4d({ 0.0, 0.0, 1.0 / std::sqrt(g_thth(r, theta)), 0.0 });
	Vector4d e3 = Vector4d({ 0.0, 0.0, 0.0, 1.0 / std::sqrt(g_pp) });

	/*frame.e1 = e1;
	frame.e2 = e2;
	frame.e3 = e3;*/

	double cy = cos(yaw), cp = cos(pitch), sy = sin(yaw), sp = sin(pitch);

	frame.e1 = (cy * cp) * e1 - (cy * sp) * e2 + sy * e3;
	frame.e2 = sp * e1 + cp * e2;
	frame.e3 = -(sy * cp) * e1 + (sp * sy) * e2 + cy * e3;

	/*frame.e1 = (sp * cy) * e1 + (sp * sy) * e2 + cp * e3;
	frame.e2 = (cp * cy) * e1 + (cp * sy) * e2 - sp * e3;
	frame.e3 = -sy * e1 + cy * e2;*/
}


void RelativisticCamera::generatePhotons() {
	double r = position[1];
	double theta = position[2];

	#pragma omp parallel for schedule(static)
	for (int i = 0; i < width * height; ++i) {
		photons.r[i] = position[1];
		photons.theta[i] = position[2];
		photons.phi[i] = position[3];
		photons.state[i] = PhotonState::Active;
		photons.dlambda[i] = 0.1;

		int x = i % width;
		int y = i / (int) width;

		// standard raytracing projection equation, first term puts (x, y) in centre of coordinates, second term applies the proper fov, x gets scaled by aspect ratio
		// must add 0.5 to offset integer screen coordinates because the constants are calculated wrong when x or y = 0 
		double screenX = (2.0 * ((double)x + 0.5) / (double)width - 1.0) * tan(fov / 2.0) * ((double)width / (double)height);
		double screenY = (1.0 - 2.0 * ((double)y + 0.5) / (double)height) * tan(fov / 2.0);

		// momenta in the camera's tetrad frame, p_t = E = 1.0
		double p_1 = 1.0 / sqrt(1.0 + screenX * screenX + screenY * screenY); // x is the forward coordinate in this instance
		double p_2 = screenY * p_1; 
		double p_3 = screenX * p_1;

		// \eta^{\mu \nu} p_\mu = -(1.0)^2 + |p_i|^2 =  -1.0 + 1.0 = 0 so it's a proper lightlike 4 vector

		// momentum projected onto global coordinates
		Vector4d p = frame.e0 + -p_1 * frame.e1 + -p_2 * frame.e2 + p_3 * frame.e3;

		// conserved covariant constants
		double E = -(p[0] * g_tt(r, theta) + p[3] * g_tphi(r, theta)); // -p_t
		double L_z = p[3] * g_phiphi(r, theta) + p[0] * g_tphi(r, theta);  // p_phi
		double Q = p[2] * p[2] * g_thth(r, theta) * g_thth(r, theta) + cos(theta) * cos(theta) * (a * a * (-E * E) + (L_z * L_z) / (sin(theta) * sin(theta)));  // carter's constant


		photons.sign_r[i] = (p[1] >= 0.0) ? 1.0f : -1.0f;
		photons.sign_theta[i] = (p[2] >= 0.0) ? 1.0f : -1.0f;
		photons.xi[i] = (L_z / E);
		photons.eta[i] = (Q / (E * E));

		PhotonDerivative derivatives = evaluate(photons.r[i], photons.theta[i], photons.sign_r[i], photons.sign_theta[i], photons.xi[i], photons.eta[i]);
		photons.dr[i] = derivatives.dr;
		photons.dtheta[i] = derivatives.dtheta;
		photons.dphi[i] = derivatives.dphi;

		/*photons.dr[i] = p[1];
		photons.dtheta[i] = p[2];
		photons.dphi[i] = p[3];*/
	}
}