#include "relativistic_camera.hpp"
#include "general_relativity.hpp"
#include "maths.hpp"

Vector4d lowerContravariant(Vector4d const& v, double r, double theta) {
	Vector4d res;
	res[0] = metricAt(0, 0, r, theta) * v[0] + metricAt(0, 3, r, theta) * v[3];
	res[1] = metricAt(1, 1, r, theta) * v[1];
	res[2] = metricAt(2, 2, r, theta) * v[2];
	res[3] = metricAt(3, 0, r, theta) * v[0] + metricAt(3, 3, r, theta) * v[3];
	return res;
}

double length2(Vector4d const& v, double r, double theta) {
	Vector4d v_lower = lowerContravariant(v, r, theta);
	double res = 0.0;
	for (size_t i = 0; i < 4; ++i)
		res += v_lower[i] * v[i];
	return res;
}

double dot(Vector4d const& v, Vector4d const& u, double r, double theta) {
	Vector4d v_lower = lowerContravariant(v, r, theta);
	double res = 0.0;
	for (size_t i = 0; i < 4; ++i)
		res += v_lower[i] * u[i];
	return res;
}

void RelativisticCamera::initTetrad() {
	//ZAMO frame 
	double r = position[1];
	double theta = position[2];

	double g_tt = metricAt(0, 0, r, theta);
	double g_rr = metricAt(1, 1, r, theta);
	double g_thth = metricAt(2, 2, r, theta);
	double g_phph = metricAt(3, 3, r, theta);
	double g_tph = metricAt(0, 3, r, theta);

	double omega = -g_tph / g_phph;
	double alpha = std::sqrt((g_tph * g_tph / g_phph) - g_tt);

	frame.e0 = Vector4d({ 1.0 / alpha, 0.0, 0.0, omega / alpha });
	Vector4d e1 = Vector4d({ 0.0, 1.0 / std::sqrt(g_rr), 0.0, 0.0 });
	Vector4d e2 = Vector4d({ 0.0, 0.0, 1.0 / std::sqrt(g_thth), 0.0 });
	Vector4d e3 = Vector4d({ 0.0, 0.0, 0.0, 1.0 / std::sqrt(g_phph) });
	double proj = g_tph * e3[3] * frame.e0[0] + g_phph * e3[3] * frame.e0[3];
	e3 = e3 + frame.e0 * proj;
	e3 /= length2(e3, r, theta);
	frame.e3 = e3;

	frame.e1 = (cos(yaw) * cos(pitch)) * e1 + (cos(yaw) * sin(pitch) * sin(roll) - sin(yaw) * cos(roll)) * e2 + (cos(yaw) * sin(pitch) * cos(roll) + sin(yaw) * sin(roll)) * e3;
	frame.e2 = (sin(yaw) * cos(pitch)) * e1 + (sin(yaw) * sin(pitch) * sin(roll) - cos(yaw) * cos(roll)) * e2 + (sin(yaw) * sin(pitch) * cos(roll) + cos(yaw) * sin(roll)) * e3;
	frame.e3 = -sin(pitch) * e1 + (cos(pitch) * sin(roll)) * e2 + (cos(pitch) * cos(roll)) * e3;
}