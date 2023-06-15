#include <cmath>
#include "../include/AzElPa.h"
#include "../include/SAT_Const.h"
#include "../include/vector.h"

//------------------------------------------------------------------------------
// void AzElPa(double s[3], double *Az, double *El, double dAds[3], double dEds[3])
//------------------------------------------------------------------------------
/**
*
* @file AzElPa.cpp
*
* Computes azimuth, elevation and partials from local tangent coordinates
*
* @param s Topocentric local tangent coordinates (East-North-Zenith frame)
* @param Az Azimuth [rad]
* @param El Elevation [rad]
* @param dAds Partials of azimuth w.r.t. s
* @param dEds Partials of elevation w.r.t. s
*/ 
//------------------------------------------------------------------------------

void AzElPa(double s[3], double *Az, double *El, double dAds[3], double dEds[3]){
	
	double pi2, rho;
	
	pi2 = 2.0*pi;
	
	rho = sqrt(s[0]*s[0]+s[1]*s[1]);
	
	// Angles
	*Az = atan2(s[0],s[1]);
	
	if (*Az<0.0){ 
		*Az = *Az+pi2;
	}
	
	*El = atan ( s[2] / rho );
	
	// Partials
	dAds[0] = s[1]/(rho*rho);
	dAds[1] = -s[0]/(rho*rho);
	dAds[2] = 0.0 ;
	
	dEds[0] = (-s[0]*s[2]/rho) / dot(s,3,s,3);
	dEds[1] = (-s[1]*s[2]/rho) / dot(s,3,s,3);
	dEds[2] = rho / dot(s,3,s,3);
}