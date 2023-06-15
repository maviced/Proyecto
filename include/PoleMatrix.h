//$Header$
//------------------------------------------------------------------------------
// PoleMatrix
//------------------------------------------------------------------------------
// EKF_GEOS3: Initial Orbit Determination using Double-R-Iteration and Extended 
// Kalman Filter methods
//
//
// Last modified:   2015/08/12   M. Mahooti
//
/**
*
* @file PoleMatrix.h
*
* Transformation from pseudo Earth-fixed to Earth-fixed coordinates for a given date
*
*/
//------------------------------------------------------------------------------

#ifndef _PoleMatrix_
#define _PoleMatrix_

void PoleMatrix(double xp, double yp, double PoleMat[3][3]);

#endif