//$Header$
//------------------------------------------------------------------------------
// NutMatrix
//------------------------------------------------------------------------------
// EKF_GEOS3: Initial Orbit Determination using Double-R-Iteration and Extended 
// Kalman Filter methods
//
//
// Last modified:   2015/08/12   M. Mahooti
//
/**
*
* @file NutMatrix.h
*
* Transformation from mean to true equator and equinox
*
*/
//------------------------------------------------------------------------------

#ifndef _NutMatrix_
#define _NutMatrix_

void NutMatrix(double Mjd_TT, double NutMat[3][3]);

#endif