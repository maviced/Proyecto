//$Header$
//------------------------------------------------------------------------------
// IERS
//------------------------------------------------------------------------------
// EKF_GEOS3: Initial Orbit Determination using Double-R-Iteration and Extended 
// Kalman Filter methods
//
//
// Last modified:   2015/08/12   M. Mahooti
//
/**
*
* @file IERS.h
*
* Management of IERS time and polar motion data
*
*/
//------------------------------------------------------------------------------

#ifndef _IERS_
#define _IERS_

void IERS(double **eop, double Mjd_TT, double *UT1_UTC, double *TAI_UTC, double *x_pole, double *y_pole);

#endif