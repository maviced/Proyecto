//$Header$
//------------------------------------------------------------------------------
// R_y
//------------------------------------------------------------------------------
// EKF_GEOS3: Initial Orbit Determination using Double-R-Iteration and Extended 
// Kalman Filter methods
//
//
// Last modified:   2015/08/12   M. Mahooti
//
/**
*
* @file R_y.h
*
* Performs a rotation around the y-axis in a three-dimensional space.
*
*/
//------------------------------------------------------------------------------

#ifndef _R_y_
#define _R_y_

void R_y(double angle, double rotmat[3][3]);

#endif