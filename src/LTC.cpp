#include <cmath>
#include "../include/LTC.h"
#include "../include/R_y.h"
#include "../include/R_z.h"
#include "../include/matriz.h"

//------------------------------------------------------------------------------
// void LTC(double lon, double lat, double M[3][3])
//------------------------------------------------------------------------------
/**
*
* @file LTC.cpp
*
* Transformation from Greenwich meridian system to local tangent coordinates
*
* @param lon Geodetic East longitude [rad]
* @param lat Geodetic latitude [rad]
* @param M Rotation matrix from the Earth equator and Greenwich meridian to the local tangent (East-North-Zenith) coordinate system
*/ 
//------------------------------------------------------------------------------


void LTC(double lon, double lat, double M[3][3]){
	
	double M1[3][3], M2[3][3], aux;
	
	R_y(-1.0*lat, M1);
	R_z(lon, M2);
	
	mult3x3(M1, M2, M);
	
	for(int i=0; i<3; i++){
		aux = M[0][i];
		M[0][i] = M[1][i];
		M[1][i] = M[2][i];
		M[2][i] = aux;
	}
}