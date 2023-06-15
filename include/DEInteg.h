//$Header$
//------------------------------------------------------------------------------
// DEInteg
//------------------------------------------------------------------------------
// EKF_GEOS3: Initial Orbit Determination using Double-R-Iteration and Extended 
// Kalman Filter methods
//
//
// Last modified:   2015/08/12   M. Mahooti
//
/**
*
* @file DEInteg.h
*
* Numerical integration methods for ordinaray differential equations.
*
* This module provides implemenation of the variable order variable 
* stepsize multistep method of Shampine & Gordon.
*
*/
//------------------------------------------------------------------------------

#ifndef _DEInteg_
#define _DEInteg_

void DEInteg(void (*func)(double, double[][1], double[][1]),double t,double tout,double relerr,double abserr,int n_eqn,double y[][1]);

#endif