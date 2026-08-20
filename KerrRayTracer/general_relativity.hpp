#pragma once
#include "maths.hpp"

// GR utility functions, since we're using carter's equations all that is necessary is the metric which is broken up so there aren't any branches, M = 1

double g_tt(double r, double theta);

double g_rr(double r, double theta);

double g_thth(double r, double theta);

double g_phiphi(double r, double theta);

double g_tphi(double r, double theta);