//$Header$
//------------------------------------------------------------------------------
// NutAngles
//------------------------------------------------------------------------------
// EKF_GEOS3: Initial Orbit Determination using Double-R-Iteration and Extended 
// Kalman Filter methods
//
//
// Last modified:   2015/08/12   M. Mahooti
//
/**
*
* @file NutAngles.h
*
* Nutation in longitude and obliquity
*
*/
//------------------------------------------------------------------------------

#ifndef _NutAngles_
#define _NutAngles_

void NutAngles(double Mjd_TT, double *dpsi, double *deps);

#endif