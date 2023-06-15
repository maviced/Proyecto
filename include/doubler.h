#include <complex>

using namespace std;

//$Header$
//------------------------------------------------------------------------------
// doubler
//------------------------------------------------------------------------------
// EKF_GEOS3: Initial Orbit Determination using Double-R-Iteration and Extended 
// Kalman Filter methods
//
//
// Last modified:   2015/08/12   M. Mahooti
//
/**
*
* @file doubler.h
*
* Performs various calculations related to double-r observations.
*
*/
//------------------------------------------------------------------------------

#ifndef _doubler_
#define _doubler_

void doubler(double cc1,double cc2,double magrsite1,double magrsite2,double magr1in,double magr2in,double los1[3],
 double los2[3],double los3[3],double rsite1[3],double rsite2[3],double rsite3[3],double t1,double t3,char direct,
 double r2[3],double r3[3], complex<double> *f1,complex<double> *f2,complex<double> *q1,double *magr1,double *magr2,double *a,double *deltae32);
 
#endif