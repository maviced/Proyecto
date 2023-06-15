//$Header$
//------------------------------------------------------------------------------
// G_AccelHarmonic
//------------------------------------------------------------------------------
// EKF_GEOS3: Initial Orbit Determination using Double-R-Iteration and Extended 
// Kalman Filter methods
//
//
// Last modified:   2015/08/12   M. Mahooti
//
/**
*
* @file G_AccelHarmonic.h
*
* Computes the gradient of the Earth's harmonic gravity field 
*
*/
//------------------------------------------------------------------------------

#ifndef _G_AccelHarmonic_
#define _G_AccelHarmonic_

void G_AccelHarmonic(double r[3][1], double U[3][3], int n_max, int m_max,double G[3][3]);

#endif