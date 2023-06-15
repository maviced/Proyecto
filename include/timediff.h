//$Header$
//------------------------------------------------------------------------------
// timediff
//------------------------------------------------------------------------------
// EKF_GEOS3: Initial Orbit Determination using Double-R-Iteration and Extended 
// Kalman Filter methods
//
//
// Last modified:   2015/08/12   M. Mahooti
//
/**
*
* @file timediff.h
*
* Time differences [s]
*
*/
//------------------------------------------------------------------------------

#ifndef _timediff_
#define _timediff_

void timediff(double UT1_UTC, double TAI_UTC, double *UT1_TAI, double *UTC_GPS, double *UT1_GPS, double *TT_UTC, double *GPS_UTC);

#endif