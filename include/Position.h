//$Header$
//------------------------------------------------------------------------------
// Position
//------------------------------------------------------------------------------
// EKF_GEOS3: Initial Orbit Determination using Double-R-Iteration and Extended 
// Kalman Filter methods
//
//
// Last modified:   2015/08/12   M. Mahooti
//
/**
*
* @file Position.h
*
* Position vector (r [m]) from geodetic coordinates (Longitude [rad], latitude [rad], altitude [m])
*
*/
//------------------------------------------------------------------------------

#ifndef _Position_
#define _Position_

void Position(double lon,double lat, double h, double r[3]);

#endif