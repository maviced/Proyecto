//$Header$
//------------------------------------------------------------------------------
// Mjday
//------------------------------------------------------------------------------
// EKF_GEOS3: Initial Orbit Determination using Double-R-Iteration and Extended 
// Kalman Filter methods
//
//
// Last modified:   2015/08/12   M. Mahooti
//
/**
*
* @file Mjday.h
*
* Modified Julian Day (MJD) calculation based on date and time components.
*
*/
//------------------------------------------------------------------------------

#ifndef _MJDAY_
#define _MJDAY_

double Mjday (int yr, int mon, int day, int hr=0, int min=0, double sec=0.0);

#endif