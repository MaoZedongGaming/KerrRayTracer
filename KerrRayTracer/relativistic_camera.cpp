#include "relativistic_camera.hpp"
#include "general_relativity.hpp"
#include "parameters.hpp"
#include "rendering_maths.hpp"
#include "maths.hpp"
#include <algorithm>
#include <iostream>
#include <omp.h>
#include <cmath>

void RelativisticCamera::turnLeft(double y) {
	yaw += y;
}

void RelativisticCamera::turnDown(double p) {
	pitch += p;
}

RelativisticCamera::RelativisticCamera(size_t w, size_t h) : width(w), height(h) {
	photons.resize(w * h);
	pixelBuffer.resize(w * h);
}

void RelativisticCamera::setPosition(double r, double theta, double phi) {
	position[0] = 0;
	position[1] = r;
	position[2] = std::clamp(theta, 0.01, PI - 0.01);
	position[3] = std::fmod(phi, TWO_PI) + (phi < 0.0) * TWO_PI;

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

	double cy = cos(yaw), cp = cos(pitch), sy = sin(yaw), sp = sin(pitch);

	frame.e1 = (cy * cp) * e1 - (cy * sp) * e2 + sy * e3;
	frame.e2 = sp * e1 + cp * e2;
	frame.e3 = -(sy * cp) * e1 + (sp * sy) * e2 + cy * e3;
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
		photons.accumulatedColour[i] = float3{ 0.0f, 0.0f, 0.0f };
		photons.transmittance[i] = 1.0f;

		int x = i % width;
		int y = i / (int) width;

		//if (x == width / 2) {
		//	photons.accumulatedColour[i] = float3{ 0.0f, 1.0f, 0.0f }; // debug pixel 
		//	photons.transmittance[i] = 0.0f;
		//}  // not the fault of the pixel at the centre of the screen

		// standard raytracing projection equation, first term puts (x, y) in centre of coordinates, second term applies the proper fov, x gets scaled by aspect ratio
		// must add 0.5 to offset integer screen coordinates because the constants are calculated wrong when x or y = 0 
		double aspectRatio = (double)width / (double)height;
		double tanHalfFov = tan(fov / 2.0);

		double screenX = (2.0 * ((double)x + 0.5) / (double)width - 1.0) * tanHalfFov * aspectRatio;
		double screenY = (1.0 - 2.0 * ((double)y + 0.5) / (double)height) * tanHalfFov;

		// momenta in the camera's tetrad frame, p_t = E = 1.0
		double p_1 = 1.0 / sqrt(1.0 + screenX * screenX + screenY * screenY); // p_1 is p_r which points OUT from the centre of the black hole
		double p_2 = screenY * p_1; 
		double p_3 = screenX * p_1;

		// \eta^{\mu \nu} p_\mu = -(1.0)^2 + |p_i|^2 =  -1.0 + 1.0 = 0 so it's a proper lightlike 4 vector

		// momentum projected onto global coordinates, REMEMBER TO FLIP p_1 BECAUSE BL COORDINATES POIMT AWAY FROM THE CENTRE!!!, flip p_2 too because e_theta points downwards
		Vector4d p = frame.e0 + -p_1 * frame.e1 + -p_2 * frame.e2 + p_3 * frame.e3;

		//if (abs(p[3]) <= 1e-14) {
		//	photons.accumulatedColour[i] = float3{ 0.0f, 1.0f, 0.0f }; // debug pixel 
		//	photons.transmittance[i] = 0.0f;
		//}

		double g_tp = g_tphi(r, theta);
		double g_pp = g_phiphi(r, theta);
		double g_t = g_tt(r, theta);
		double g_tt = g_thth(r, theta);
		double sinTh = sin(theta);
		double cosTh = cos(theta);

		// conserved covariant constants
		double E = -(p[0] * g_t + p[3] * g_tp); // -p_t
		double L_z = p[3] * g_pp + p[0] * g_tp;  // p_phi
		double Q = p[2] * p[2] * g_tt * g_tt + cosTh * cosTh * (a * a * (-E * E) + (L_z * L_z) / (sinTh * sinTh));  // carter's constant


		photons.sign_r[i] = (p[1] >= 0.0) ? 1.0f : -1.0f;
		photons.sign_theta[i] = (p[2] >= 0.0) ? 1.0f : -1.0f;
		photons.xi[i] = (L_z / E);
		photons.eta[i] = (Q / (E * E));

		PhotonDerivative derivatives = evaluate(photons.r[i], photons.theta[i], photons.sign_r[i], photons.sign_theta[i], photons.xi[i], photons.eta[i]);
		photons.dr[i] = derivatives.dr;
		photons.dtheta[i] = derivatives.dtheta;
		photons.dphi[i] = derivatives.dphi;

		// don't add noise because the the conserved constants are finicky 
		/*photons.dr[i] = p[1];
		photons.dtheta[i] = p[2];
		photons.dphi[i] = p[3];*/
	}
}