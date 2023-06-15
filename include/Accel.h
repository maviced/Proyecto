//$Header$
//------------------------------------------------------------------------------
// Accel
//------------------------------------------------------------------------------
// EKF_GEOS3: Initial Orbit Determination using Double-R-Iteration and Extended 
// Kalman Filter methods
//
//
// Last modified:   2015/08/12   M. Mahooti
//
/**
*
* @file Accel.h
*
*   Computes the acceleration of an Earth orbiting satellite due to 
*    - the Earth's harmonic gravity field, 
*    - the gravitational perturbations of the Sun and Moon
*    - the solar radiation pressure and
*    - the atmospheric drag
*
*/
//------------------------------------------------------------------------------

#ifndef _Accel_
#define _Accel_

void Accel(double x,double Y[6][1], double dY[6][1]);

#endif