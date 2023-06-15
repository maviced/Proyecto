//$Header$
//------------------------------------------------------------------------------
// MeasUpdate
//------------------------------------------------------------------------------
// EKF_GEOS3: Initial Orbit Determination using Double-R-Iteration and Extended 
// Kalman Filter methods
//
//
// Last modified:   2015/08/12   M. Mahooti
//
/**
*
* @file MeasUpdate.h
*
* Implements a measurement update step for a Kalman filter
*
*/
//------------------------------------------------------------------------------

#ifndef _MeasUpdate_
#define _MeasUpdate_

void MeasUpdate(double x[6][1], double  z, double g, double s,double G[1][6], double P[6][6],int n, double K[6][1], double x2[6][1],double P2[6][6]);

#endif