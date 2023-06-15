#ifndef _Legendre_
#define _Legendre_

//$Header$
//------------------------------------------------------------------------------
// Legendre
//------------------------------------------------------------------------------
// EKF_GEOS3: Initial Orbit Determination using Double-R-Iteration and Extended 
// Kalman Filter methods
//
//
// Last modified:   2015/08/12   M. Mahooti
//
/**
*
* @file Legendre.h
*
* Implements the computation of Legendre polynomials and their derivatives using a recurrence relation.
*
*/
//------------------------------------------------------------------------------

void Legendre(int n,int m,double fi, double pnm[362][362], double dpnm[362][362]);

#endif