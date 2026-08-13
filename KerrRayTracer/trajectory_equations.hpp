//no need to separate lightlike and timelike geodesics, since the equations are the same, just with different initial conditions

#include "maths.hpp"
#include <cmath>
#include "length_scales.hpp"

Vector4d timelikeTrajectoryEquations(Vector4d const& pos, Vector4d const& momentum, float particleMass, float blackHoleMass, float spin) {
	double r = pos[1];
	double theta = pos[2];
	double p_r = momentum[1];
	double p_theta = momentum[2];
	double angularMomentum = momentum[3];
	double energy = momentum[0];
	
	double a = spinParameter(spin, blackHoleMass);
	double deltaVal = delta(r, blackHoleMass, spin);
	double sigmaVal = sigma(r, theta, blackHoleMass, spin);
	double mu2 = particleMass * particleMass;
	double carterConstant = p_theta * p_theta + cos(theta) * cos(theta) * (a * a * (mu2 - energy * energy) + (angularMomentum * angularMomentum) / (sin(theta) * sin(theta)));

	double P = energy * (r * r + a * a) - a * angularMomentum;
	double R = P * P - deltaVal * (mu2 * r * r + (angularMomentum - a * energy) * (angularMomentum - a * energy) + carterConstant);
	double Theta = carterConstant - cos(theta) * cos(theta) * (a * a * (mu2 - energy * energy) + (angularMomentum * angularMomentum) / (sin(theta) * sin(theta)));

	Vector4d derivatives;
	
	derivatives[0] = (a * (angularMomentum - a * energy * sin(theta) * sin(theta)) + (r * r + a * a) * P / deltaVal) / sigmaVal;
	derivatives[1] = sqrt(R) / sigmaVal;
	derivatives[2] = std::abs(sqrt(Theta)) / sigmaVal;
	derivatives[3] = (angularMomentum/(sin(theta) * sin(theta)) - a * energy + a * P / deltaVal) / sigmaVal;
	return derivatives;
}

Vector4d lightlikeTrajectoryEquations(Vector4d const& pos, Vector4d const& momentum, float particleMass, float blackHoleMass, float spin) {
	double r = pos[1];
	double theta = pos[2];
	double p_r = momentum[1];
	double p_theta = momentum[2];
	double angularMomentum = momentum[3];
	double energy = momentum[0];

	double a = spinParameter(spin, blackHoleMass);
	double deltaVal = delta(r, blackHoleMass, spin);
	double sigmaVal = sigma(r, theta, blackHoleMass, spin);
	double carterConstant = p_theta * p_theta + cos(theta) * cos(theta) * (a * a * (- energy * energy) + (angularMomentum * angularMomentum) / (sin(theta) * sin(theta)));

	double P = energy * (r * r + a * a) - a * angularMomentum;
	double R = P * P - deltaVal * ((angularMomentum - a * energy) * (angularMomentum - a * energy) + carterConstant);
	double Theta = carterConstant - cos(theta) * cos(theta) * (a * a * (- energy * energy) + (angularMomentum * angularMomentum) / (sin(theta) * sin(theta)));

	Vector4d derivatives;

	derivatives[0] = (a * (angularMomentum - a * energy * sin(theta) * sin(theta)) + (r * r + a * a) * P / deltaVal) / sigmaVal;
	derivatives[1] = sqrt(R) / sigmaVal;
	derivatives[2] = std::abs(sqrt(Theta)) / sigmaVal;
	derivatives[3] = (angularMomentum / (sin(theta) * sin(theta)) - a * energy + a * P / deltaVal) / sigmaVal;
	return derivatives;
}