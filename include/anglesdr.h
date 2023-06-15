//$Header$
//------------------------------------------------------------------------------
// anglesdr
//------------------------------------------------------------------------------
// EKF_GEOS3: Initial Orbit Determination using Double-R-Iteration and Extended 
// Kalman Filter methods
//
//
// Last modified:   2015/08/12   M. Mahooti
//
/**
*
* @file anglesdr.h
*
* This function solves the problem of orbit determination using three optical sightings.
*
*/
//------------------------------------------------------------------------------

#ifndef _anglesdr_
#define _anglesdr_

void anglesdr ( double az1, double az2,double az3, double el1, double el2, double el3, double Mjd1, double Mjd2, double Mjd3, double rsite1[3], 
 double rsite2[3], double rsite3[3], double r2[3], double v2[3] );

#endif