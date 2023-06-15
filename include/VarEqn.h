//$Header$
//------------------------------------------------------------------------------
// VarEqn
//------------------------------------------------------------------------------
// EKF_GEOS3: Initial Orbit Determination using Double-R-Iteration and Extended 
// Kalman Filter methods
//
//
// Last modified:   2015/08/12   M. Mahooti
//
/**
*
* @file VarEqn.h
*
* Computes the variational equations, i.e. the derivative of the state vector
*   and the state transition matrix
*
*/
//------------------------------------------------------------------------------

#ifndef _VarEqn_
#define _VarEqn_

void VarEqn(double x, double yPhi[42][1], double yPhip[42][1]);

#endif