#pragma once

constexpr bool ENABLE_OUTER_HORIZON = true;
constexpr bool ENABLE_ACCRETION_DISK = true;
constexpr bool ENABLE_DOPPLER_BEAMING = true;
constexpr bool ENABLE_DISK_DENSITY = true;
constexpr bool ENABLE_OPAQUE_DISK = false;
constexpr bool PIXEL_DEBUG = false;
constexpr bool RKDP_INTEGRATION = false;  // 30k steps is good
constexpr bool RK4_INTEGRATION = true;  // fuuuuck heuristic RK4 is just as good, just use heuristic RK4 
constexpr bool RK5_INTEGRATION = false; // not any better than RK4