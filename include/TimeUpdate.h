//$Header$
//------------------------------------------------------------------------------
// TimeUpdate
//------------------------------------------------------------------------------
// EKF_GEOS3: Initial Orbit Determination using Double-R-Iteration and Extended 
// Kalman Filter methods
//
//
// Last modified:   2015/08/12   M. Mahooti
//
/**
*
* @file TimeUpdate.h
*
* Implements a time update step for a Kalman filter.
*
*/
//------------------------------------------------------------------------------

#ifndef _TimeUpdate_
#define _TimeUpdate_

void TimeUpdate(double P[6][6], double Phi[6][6], double Qdt=0.0);

#endif