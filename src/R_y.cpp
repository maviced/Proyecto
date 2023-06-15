#include <cmath>
#include "../include/R_y.h"

//------------------------------------------------------------------------------
// void R_y(double angle, double rotmat[3][3])
//------------------------------------------------------------------------------
/**
*
* @file R_y.cpp
*
* Performs a rotation around the y-axis in a three-dimensional space.
*
* @param angle angle of rotation [rad]
* @param rotmat vector result
*/ 
//------------------------------------------------------------------------------

void R_y(double angle, double rotmat[3][3]) {
	
	double C, S;
    C = cos(angle);
    S = sin(angle);
    rotmat[0][0] = C; rotmat[0][1] = 0.0; rotmat[0][2] = -1.0*S;
    rotmat[1][0] = 0.0; rotmat[1][1] = 1.0; rotmat[1][2] = 0.0;
    rotmat[2][0] = S; rotmat[2][1] = 0.0; rotmat[2][2] = C;
}