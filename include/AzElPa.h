//$Header$
//------------------------------------------------------------------------------
// AzElPa
//------------------------------------------------------------------------------
// EKF_GEOS3: Initial Orbit Determination using Double-R-Iteration and Extended 
// Kalman Filter methods
//
//
// Last modified:   2015/08/12   M. Mahooti
//
/**
*
* @file AzElPa.h
*
* Computes azimuth, elevation and partials from local tangent coordinates
*
*/
//------------------------------------------------------------------------------

#ifndef _AzEIPa_
#define _AzEIPa_

void AzElPa(double s[3], double *Az, double *El, double dAds[3], double dEds[3]);

#endif