//$Header$
//------------------------------------------------------------------------------
// AccelHarmonic
//------------------------------------------------------------------------------
// EKF_GEOS3: Initial Orbit Determination using Double-R-Iteration and Extended 
// Kalman Filter methods
//
//
// Last modified:   2015/08/12   M. Mahooti
//
/**
*
* @file AccelHarmonic.h
*
* Computes the acceleration due to the harmonic gravity field of the central body
*
*/
//------------------------------------------------------------------------------

#ifndef _AccelHarmonic_
#define _AccelHarmonic_

void AccelHarmonic(double r[3][1], double E[3][3], int n_max, int m_max,double a[3][1]);

#endif