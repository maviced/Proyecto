#include <cmath>
#include <stdexcept>
#include "../include/Geodetic.h"
#include "../include/SAT_Const.h"
#include "../include/vector.h"

//------------------------------------------------------------------------------
// void Geodetic( double r[3], double *lon, double *lat, double *h)
//------------------------------------------------------------------------------
/**
*
* @file Geodetic.cpp
*
* geodetic coordinates (Longitude [rad], latitude [rad], altitude [m]) from given position vector (r [m])
*
* @param r given position vector (r [m])
* @param lon Longitude [rad] from vector r
* @param lat latitude [rad] from vector r
* @param h altitude [m] from vector r
*/ 
//------------------------------------------------------------------------------

void Geodetic( double r[3], double *lon, double *lat, double *h){
	
	double R_equ, f, epsRequ, e2, X, Y, Z, rho2, dZ, ZdZ, Nh, SinPhi, N, dZ_new;
	
	R_equ = R_Earth;
	f     = f_Earth;

	epsRequ = 2.220446049250313e-16*R_equ;        // Convergence criterion
	e2      = f*(2.0-f);        // Square of eccentricity
	
	X = r[0];                   // Cartesian coordinates
	Y = r[1];
	Z = r[2];
	rho2 = X*X + Y*Y;           // Square of distance from z-axis
	
	// Check validity of input data
	if (norm(r, 3)==0.0){
		throw std::runtime_error("Invalid input in Geodetic constructor");
		*lon = 0.0;
		*lat = 0.0;
		*h   = -R_Earth;
	}
	
	// Iteration 
	dZ = e2*Z;
	
	while(true){
		ZdZ    =  Z + dZ;
		Nh     =  sqrt ( rho2 + ZdZ*ZdZ ); 
		SinPhi =  ZdZ / Nh;                    // Sine of geodetic latitude
		N      =  R_equ / sqrt(1.0-e2*SinPhi*SinPhi);
		dZ_new =  N*e2*SinPhi;
		if ( fabs(dZ-dZ_new) < epsRequ ){
			break;
		}
		dZ = dZ_new;
	}
	
	// Longitude, latitude, altitude
	*lon = atan2 ( Y, X );
	*lat = atan2 ( ZdZ, sqrt(rho2) );
	*h   = Nh - N;

}