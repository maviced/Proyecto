//------------------------------------------------------------------------------
// void PoleMatrix(double xp, double yp, double PoleMat[3][3])
//------------------------------------------------------------------------------
/**
*
* @file PoleMatrix.cpp
*
* Transformation from pseudo Earth-fixed to Earth-fixed coordinates for a given date
*
* @param xp Pole coordinte
* @param yp Pole coordinte
* @param PoleMat Pole matrix
*/ 
//------------------------------------------------------------------------------

#include <cmath>
#include "../include/PoleMatrix.h"
#include "../include/R_y.h"
#include "../include/R_x.h"
#include "../include/matriz.h"

void PoleMatrix(double xp, double yp, double PoleMat[3][3]){
	
	double m1[3][3], m2[3][3];
	
	R_y(-xp,m1);
	R_x(-yp, m2);
	
	mult3x3(m1, m2, PoleMat);
}