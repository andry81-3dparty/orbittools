//
// globals.h
//
// Copyright (c) 2003-2013 Michael F. Henry
//
#pragma once

#include "orbitTools/Defines.h"

#include "orbitTools/SysUndefine.h"

#include <cmath>
#include <cfloat>
#include <limits> // can't use DBL_MAX in gcc, won't compile because of an error: `error: call of overloaded 'dd_real(long double)' is ambiguous`

#include "orbitTools/SysDefine.h"


namespace Zeptomoby 
{
namespace OrbitTools
{
#if QD_INTEGRATION_ENABLED
const double PI           = double::_pi(); //3.141592653589793;
const double TWOPI        = double::_2pi(); //2.0 * PI;
#else
const double PI           = 3.141592653589793;
const double TWOPI        = 2.0 * PI;
#endif
const double RADS_PER_DEG = PI / 180.0;

const double GM           = double(398601'2) / 10;       // Earth gravitational constant, km^3/sec^2
const double GEOSYNC_ALT  = double(42241'892) / 1000;    // km
const double EARTH_DIA    = double(12742'0176) / 10000;  //12800.0;      // km
const double DAY_SIDERAL  = (23 * 3600) + (56 * 60) + double(409) / 100;  // sec
const double DAY_24HR     = (24 * 3600);   // sec

const double AE           = 1.0;
const double AU           = double(149597870'691) / 1000; //149597870.0;  // Astronomical unit (km) (IAU 76)
const double SR           = 696000.0;     // Solar radius (km)      (IAU 76)
const double XKMPER_WGS72 = double(6378'137) / 1000; //6378.135;            // Earth equatorial radius - km (WGS '72)
const double F            = (double(6378'1370) - 6356'7523) / 6378'1370; //1.0 / 298.26; // Earth flattening (WGS '72)
const double GE           = double(398600'8) / 10;      // Earth gravitational constant (WGS '72)
const double J2           = double(1'0826158) / 1e10;    // J2 harmonic (WGS '72)
const double J3           = double(-2'53881) / 1e11;     // J3 harmonic (WGS '72)
const double J4           = double(-1'65597) / 1e11;     // J4 harmonic (WGS '72)
const double CK2          = J2 / 2.0;
const double CK4          = -3.0 * J4 / 8.0;
const double XJ3          = J3;
const double QO           = AE + 120.0 / XKMPER_WGS72;
const double S            = AE + 78.0  / XKMPER_WGS72;
const double HR_PER_DAY   = 24.0;          // Hours per day   (solar)
const double MIN_PER_DAY  = 1440.0;        // Minutes per day (solar)
const double SEC_PER_DAY  = 86400.0;       // Seconds per day (solar)
const double OMEGA_E      = double(1'00273790934) / 1e11; // earth rotation per sideral day
const double XKE          = std::sqrt(3600.0 * GE /           //sqrt(ge) ER^3/min^2
                                (XKMPER_WGS72 * XKMPER_WGS72 * XKMPER_WGS72)); 
const double QOMS2T       = std::pow((QO - S), 4);            //(QO - S)^4 ER^4

// `double_type` instead `double` to bypass definition substitution
const double_type double_max = (std::numeric_limits<double_type>::max)();

// Utility functions
double sqr   (const double x);
double Fmod2p(const double arg);
double AcTan (const double sinx, const double cosx);

double rad2deg(const double);
double deg2rad(const double);

// TO FIX triginometric range in acos/asin functions !!!

extern inline double truncate_float_to_minmax(double value, double min_value, double max_value)
{
    if (max_value < value) {
        return max_value;
    }

    if (value < min_value) {
        return min_value;
    }

    return value;
}

extern inline double fix_float_trigonometric_range_factor(double value)
{
    // avoid fix in special case
    if (std::isnormal(value) && value != double_max && value != -double_max) {
        return truncate_float_to_minmax(value, -1.0, +1.0);
    }

    return value;
}

}
}

#include "orbitTools/SysUndefine.h"
