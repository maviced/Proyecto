#include <cmath>
#include "../include/GHAMatrix.h"
#include "../include/R_z.h"
#include "../include/gast.h"

//------------------------------------------------------------------------------
// void GHAMatrix(double Mjd_UT1, double GHAmat[3][3])
//------------------------------------------------------------------------------
/**
*
* @file GHAMatrix.cpp
*
* Transformation from true equator and equinox to Earth equator and Greenwich meridian system 
*
* @param Mjd_UT1 Modified Julian Date UT1
* @return Greenwich Hour Angle matrix
*/ 
//------------------------------------------------------------------------------

void GHAMatrix(double Mjd_UT1, double GHAmat[3][3]){
	
	R_z(gast(Mjd_UT1), GHAmat);
}