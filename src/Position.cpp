#include <cmath>
#include "../include/Position.h"
#include "../include/SAT_Const.h"

//------------------------------------------------------------------------------
// void Position(double lon,double lat, double h, double r[3])
//------------------------------------------------------------------------------
/**
*
* @file Position.cpp
*
* Position vector (r [m]) from geodetic coordinates (Longitude [rad], latitude [rad], altitude [m])
*
* @param lon Longitude [rad
* @param lat latitude [rad]
* @param h altitude [m]
* @param r Position vector (r [m]) from geodetic coordinates
*/ 
//------------------------------------------------------------------------------

void Position(double lon,double lat, double h, double r[3]){
	
	double R_equ, f, e2, CosLat, SinLat, N;
	
	R_equ = R_Earth;
	f = f_Earth;
	
	e2 = f*(2.0-f);      // Square of eccentricity
	CosLat = cos(lat);   // (Co)sine of geodetic latitude
	SinLat = sin(lat);
	
	// Position vector
	N = R_equ / sqrt(1.0-e2*SinLat*SinLat);
	
	r[0] =  (         N+h)*CosLat*cos(lon); 
	r[1] =  (         N+h)*CosLat*sin(lon);
	r[2] =  ((1.0-e2)*N+h)*SinLat;
}