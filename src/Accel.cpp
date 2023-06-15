#include <cmath>
#include <iostream>
#include "../include/Accel.h"
#include "../include/AccelHarmonic.h"
#include "../include/SAT_Const.h"
#include "../include/IERS.h"
#include "../include/timediff.h"
#include "../include/globales.h"
#include "../include/PrecMatrix.h"
#include "../include/NutMatrix.h"
#include "../include/matriz.h"
#include "../include/PoleMatrix.h"
#include "../include/GHAMatrix.h"

using namespace std;

//------------------------------------------------------------------------------
// void Accel(double x,double Y[6][1], double dY[6][1])
//------------------------------------------------------------------------------
/**
*
* @file Accel.cpp
*
*   Computes the acceleration of an Earth orbiting satellite due to 
*    - the Earth's harmonic gravity field, 
*    - the gravitational perturbations of the Sun and Moon
*    - the solar radiation pressure and
*    - the atmospheric drag
*
* @param x Terrestrial Time (Modified Julian Date)
* @param Y Satellite state vector in the ICRF/EME2000 system
* @param dY Acceleration (a=d^2r/dt^2) in the ICRF/EME2000 system
*/ 
//------------------------------------------------------------------------------

void Accel(double x,double Y[6][1], double dY[6][1]){
	
	double UT1_UTC, TAI_UTC, x_pole, y_pole, UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC, Mjd_UT1, P[3][3], N[3][3], T[3][3], E1[3][3], E2[3][3], E[3][3],
	a[3][1];
	
	IERS(eopdata, AuxParam.Mjd_TT + x/86400, &UT1_UTC, &TAI_UTC, &x_pole, &y_pole);
	timediff(UT1_UTC, TAI_UTC, &UT1_TAI, &UTC_GPS, &UT1_GPS, &TT_UTC, &GPS_UTC);
	
	Mjd_UT1 = AuxParam.Mjd_TT + x/86400 + (UT1_UTC-TT_UTC)/86400.0;
	
	PrecMatrix(MJD_J2000,AuxParam.Mjd_TT + x/86400, P);
	NutMatrix(AuxParam.Mjd_TT + x/86400, N);
	mult3x3(N, P, T);
	PoleMatrix(x_pole,y_pole, E1);
	GHAMatrix(Mjd_UT1, E2);
	mult3Matrices(E1, E2, T, E);
	
	double Y2[3][1] = {{Y[0][0]}, {Y[1][0]}, {Y[2][0]}};
	
	AccelHarmonic(Y2, E, AuxParam.n, AuxParam.m, a);
	
	dY[0][0] = Y[3][0];
	dY[1][0] = Y[4][0];
	dY[2][0] = Y[5][0];
	dY[3][0] = a[0][0];
	dY[4][0] = a[1][0];
	dY[5][0] = a[2][0];
	
}