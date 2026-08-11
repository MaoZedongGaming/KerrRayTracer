#pragma once
#include "maths.hpp"
#include "length_scales.hpp"

//metric stuff, christoffel symbols, etc. for Kerr spacetime in Boyer-Lindquist coordinates

template <typename T>
T metricAt(size_t mu, size_t nu, Vector4d const& pos, T mass, T spin);

double metricDerivativeAt(size_t alpha, size_t mu, size_t nu, Vector4d const& pos, double mass, double spin);

double inverseMetricAt(size_t mu, size_t nu, Vector4d const& pos, double mass, double spin);

double christoffelSymbol(size_t lambda, size_t mu, size_t nu, Vector4d const& pos, double mass, double spin);