//$Header$
//------------------------------------------------------------------------------
// PrecMatrix
//------------------------------------------------------------------------------
// EKF_GEOS3: Initial Orbit Determination using Double-R-Iteration and Extended 
// Kalman Filter methods
//
//
// Last modified:   2015/08/12   M. Mahooti
//
/**
*
* @file PrecMatrix.h
*
* Precession transformation of equatorial coordinates
*
*/
//------------------------------------------------------------------------------

#ifndef _PrecMatrix_
#define _PrecMatrix_

void PrecMatrix(double Mjd_1, double Mjd_2, double PrecMat[3][3]);

#endif