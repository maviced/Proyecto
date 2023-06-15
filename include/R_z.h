//$Header$
//------------------------------------------------------------------------------
// R_z
//------------------------------------------------------------------------------
// EKF_GEOS3: Initial Orbit Determination using Double-R-Iteration and Extended 
// Kalman Filter methods
//
//
// Last modified:   2015/08/12   M. Mahooti
//
/**
*
* @file R_z.h
*
* Performs a rotation around the z-axis in a three-dimensional space.
*
*/
//------------------------------------------------------------------------------

#ifndef _R_z_
#define _R_z_

void R_z(double angle, double rotmat[3][3]);

#endif