#pragma once

constexpr bool ENABLE_ACCRETION_DISK = true;
constexpr bool ENABLE_DOPPLER_BEAMING = true;
constexpr bool PIXEL_DEBUG = false;
constexpr bool RKDP_INTEGRATION = true;  // 30k steps is good
constexpr bool RK4_INTEGRATION = false;  // fuuuuck adaptive RK4 is just as good
constexpr bool RK5_INTEGRATION = false; // not any better than RK4