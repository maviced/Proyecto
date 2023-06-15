#include <cmath>
#include "../include/NutMatrix.h"
#include "../include/MeanObliquity.h"
#include "../include/NutAngles.h"
#include "../include/R_x.h"
#include "../include/R_z.h"
#include "../include/matriz.h"

//------------------------------------------------------------------------------
// void NutMatrix(double Mjd_TT, double NutMat[3][3])
//------------------------------------------------------------------------------
/**
*
* @file NutMatrix.cpp
*
* Transformation from mean to true equator and equinox
*
* @param Mjd_TT Modified Julian Date (Terrestrial Time)
* @param NutMat Nutation matrix
*/ 
//------------------------------------------------------------------------------

void NutMatrix(double Mjd_TT, double NutMat[3][3]){
	
	double eps, dpsi, deps, r1[3][3], r2[3][3], r3[3][3];
	
	// Mean obliquity of the ecliptic
	eps = MeanObliquity (Mjd_TT);

	// Nutation in longitude and obliquity
	NutAngles (Mjd_TT, &dpsi, &deps);

	// Transformation from mean to true equator and equinox
	R_x(-eps-deps, r1);
	R_z(-dpsi, r2);
	R_x(+eps, r3);
	mult3Matrices(r1,r2,r3,NutMat);
}