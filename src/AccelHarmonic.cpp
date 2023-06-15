#include <cmath>
#include <iostream>
#include "../include/AccelHarmonic.h"
#include "../include/SAT_Const.h"
#include "../include/matriz.h"
#include "../include/vector.h"
#include "../include/Legendre.h"

//------------------------------------------------------------------------------
// void AccelHarmonic(double r[3][1], double E[3][3], int n_max, int m_max,double a[3][1])
//------------------------------------------------------------------------------
/**
*
* @file AccelHarmonic.cpp
*
* Computes the acceleration due to the harmonic gravity field of the central body
*
* @param Mjd_TT Modified Julian Date of TT
* @param r Satellite position vector in the inertial system
* @param E Transformation matrix to body-fixed system
* @param n_max Maximum degree
* @param m_max Maximum order (m_max<=n_max; m_max=0 for zonals, only)
* @param a Acceleration (a=d^2r/dt^2)
*/ 
//------------------------------------------------------------------------------

using namespace std;

extern double Cnm[362][362], Snm[362][362], pnm[362][362], dpnm[362][362];

void AccelHarmonic(double r[3][1], double E[3][3], int n_max, int m_max,double a[3][1]){
	
	double gm, r_ref, r_bf[3][1], v1[3], d, latgc, lon, dUdr, dUdlatgc, dUdlon, q3, q2, q1, b1, b2, b3, r2xy, ax, ay, az, a_bf[3][1];
	
	gm    = 398600.4415e9;              // [m^3/s^2]; JGM3/EGM96
	r_ref = 6378.1363e3;                // Radius Earth [m]; JGM3/EGM96
	
	// Body-fixed position 
	mult3x1(E, r, r_bf);
	
	// Auxiliary quantities
	matrizAVector(r_bf, v1);
	d = norm(v1, 3);                     // distance
	latgc = asin(r_bf[2][0]/d);
	lon = atan2(r_bf[1][0],r_bf[0][0]);
	
	Legendre(n_max,m_max,latgc, pnm, dpnm);
	
	dUdr = 0;
	dUdlatgc = 0;
	dUdlon = 0;
	q3 = 0; q2 = q3; q1 = q2;
	
	
	for (int n=0; n<n_max+1; n++){
		
		b1 = (-gm/(d*d))*pow(r_ref/d, n)*(n+1);
		b2 =  (gm/d)*pow(r_ref/d,n);
		b3 =  (gm/d)*pow(r_ref/d,n);
		
		for (int m=0; m<n+1; m++){
			q1 = q1 + pnm[n][m]*(Cnm[n][m]*cos(m*lon)+Snm[n][m]*sin(m*lon));
			q2 = q2 + dpnm[n][m]*(Cnm[n][m]*cos(m*lon)+Snm[n][m]*sin(m*lon));
			q3 = q3 + m*pnm[n][m]*(Snm[n][m]*cos(m*lon)-Cnm[n][m]*sin(m*lon));
		}
		dUdr     = dUdr     + q1*b1;
		dUdlatgc = dUdlatgc + q2*b2;
		dUdlon   = dUdlon   + q3*b3;
		q3 = 0; q2 = q3; q1 = q2;
	}
	
	// Body-fixed acceleration
	r2xy = r_bf[0][0]*r_bf[0][0]+r_bf[1][0]*r_bf[1][0];
	
	ax = (1/d*dUdr-r_bf[2][0]/(d*d*sqrt(r2xy))*dUdlatgc)*r_bf[0][0]-(1/r2xy*dUdlon)*r_bf[1][0];
	ay = (1/d*dUdr-r_bf[2][0]/(d*d*sqrt(r2xy))*dUdlatgc)*r_bf[1][0]+(1/r2xy*dUdlon)*r_bf[0][0];
	az =  1/d*dUdr*r_bf[2][0]+sqrt(r2xy)/(d*d)*dUdlatgc;
	
	a_bf[0][0] = ax;
	a_bf[1][0] = ay;
	a_bf[2][0] = az;
	
	// Inertial acceleration 
	transpuesta(E, 3, 3);
	mult3x1(E, a_bf, a);
}