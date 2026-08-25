#include "general_relativity.hpp"
#include "parameters.hpp"
//
//double g_tt(double r, double theta) {
//	return -(1 - 2 * r / sigma(r, theta));
//}
//
//double g_rr(double r, double theta) {
//	return sigma(r, theta)/delta(r);
//}
//
//double g_thth(double r, double theta) {
//	return sigma(r, theta);
//}
//
//double g_phiphi(double r, double theta) {
//	return (r * r + a * a + 2 * r * a * a * std::sin(theta) * std::sin(theta) / sigma(r, theta)) * std::sin(theta) * std::sin(theta);
//}
//
//double g_tphi(double r, double theta) {
//	return -4 * r * a * std::sin(theta) * std::sin(theta) / sigma(r, theta);
//}


double g_tt(double r, double theta) {
    double r2 = r * r;
    double cosTh = cos(theta);
    double sigma = r2 + a * a * cosTh * cosTh;
    return -(1.0 - (2.0 * r) / sigma);
}

double g_rr(double r, double theta) {
    double r2 = r * r;
    double a2 = a * a;
    double cosTh = cos(theta);
    double sigma = r2 + a2 * cosTh * cosTh;
    double delta = r2 - 2 * r + a2;
    return sigma / delta;
}

double g_tphi(double r, double theta) {
    double r2 = r * r;
    double sinTh = sin(theta);
    double cosTh = cos(theta);
    double sigma = r2 + a * a * cosTh * cosTh;
    return -(4.0 * a * r * sinTh * sinTh) / sigma;
}

double g_phiphi(double r, double theta) {
    double r2 = r * r;
    double a2 = a * a;
    double sinTh = sin(theta);
    double sin2 = sinTh * sinTh;
    double cosTh = cos(theta);
    double sigma = r2 + a2 * cosTh * cosTh;
    return (r2 + a2 + (2.0 * a2 * r * sin2) / sigma) * sin2;
}

double g_thth(double r, double theta) {
    double r2 = r * r;
    double cosTh = cos(theta);
    return r2 + a * a * cosTh * cosTh;
}
