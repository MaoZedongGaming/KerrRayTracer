#pragma once
#include "maths.hpp"

//metric stuff, christoffel symbols, etc. for Kerr spacetime in Boyer-Lindquist coordinates
//kinda a useless file because carter's equations are overpowered and don't require the metric or christoffel symbols, but it's here for completeness
// still need the metric and inverse metric for tetrad and gravitational redshift calculations and raising/lowering indices, so it's not completely useless

double metricAt(size_t mu, size_t nu, double r, double theta);  // kerr metric doesn't depend on t or phi, so we can just pass in r and theta instead of the full position vector, also M = 1 

//double metricDerivativeAt(size_t alpha, size_t mu, size_t nu, Vector4d const& pos, double mass, double spin);

double inverseMetricAt(size_t mu, size_t nu, double r, double theta);

//double christoffelSymbol(size_t lambda, size_t mu, size_t nu, Vector4d const& pos, double mass, double spin);