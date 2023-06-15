//$Header$
//------------------------------------------------------------------------------
// LTC
//------------------------------------------------------------------------------
// EKF_GEOS3: Initial Orbit Determination using Double-R-Iteration and Extended 
// Kalman Filter methods
//
//
// Last modified:   2015/08/12   M. Mahooti
//
/**
*
* @file LTC.h
*
* Transformation from Greenwich meridian system to local tangent coordinates
*
*/
//------------------------------------------------------------------------------

#ifndef _LTC_
#define _LTC_

void LTC(double lon, double lat, double M[3][3]);

#endif