//$Header$
//------------------------------------------------------------------------------
// R_x
//------------------------------------------------------------------------------
// EKF_GEOS3: Initial Orbit Determination using Double-R-Iteration and Extended 
// Kalman Filter methods
//
//
// Last modified:   2015/08/12   M. Mahooti
//
/**
*
* @file R_x.h
*
* Performs a rotation around the x-axis in a three-dimensional space.
*
*/
//------------------------------------------------------------------------------

#ifndef _R_x_
#define _R_x_

void R_x(double angle, double rotmat[3][3]);

#endif