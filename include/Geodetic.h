//$Header$
//------------------------------------------------------------------------------
// Geodetic
//------------------------------------------------------------------------------
// EKF_GEOS3: Initial Orbit Determination using Double-R-Iteration and Extended 
// Kalman Filter methods
//
//
// Last modified:   2015/08/12   M. Mahooti
//
/**
*
* @file Geodetic.h
*
* geodetic coordinates (Longitude [rad], latitude [rad], altitude [m]) from given position vector (r [m])
*
*/
//------------------------------------------------------------------------------

#ifndef _Geodetic_
#define _Geodetic_

void Geodetic( double r[3], double *lon, double *lat, double *h);

#endif