//$Header$
//------------------------------------------------------------------------------
// GHAMatrix
//------------------------------------------------------------------------------
// EKF_GEOS3: Initial Orbit Determination using Double-R-Iteration and Extended 
// Kalman Filter methods
//
//
// Last modified:   2015/08/12   M. Mahooti
//
/**
*
* @file GHAMatrix.h
*
* Transformation from true equator and equinox to Earth equator and Greenwich meridian system 
*
*/
//------------------------------------------------------------------------------

#ifndef _GHAMatrix_
#define _GHAMatrix_

void GHAMatrix(double Mjd_UT1, double GHAmat[3][3]);

#endif