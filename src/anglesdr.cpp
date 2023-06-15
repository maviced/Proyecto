#include <cmath>
#include <complex>
#include <iostream>
#include <iomanip>
#include "../include/anglesdr.h"
#include "../include/SAT_Const.h"
#include "../include/Geodetic.h"
#include "../include/LTC.h"
#include "../include/matriz.h"
#include "../include/IERS.h"
#include "../include/timediff.h"
#include "../include/NutMatrix.h"
#include "../include/PrecMatrix.h"
#include "../include/PoleMatrix.h"
#include "../include/GHAMatrix.h"
#include "../include/vector.h"
#include "../include/doubler.h"

//------------------------------------------------------------------------------
// void anglesdr ( double az1, double az2,double az3, double el1, double el2, double el3, double Mjd1, double Mjd2, double Mjd3, double rsite1[3], 
//	double rsite2[3], double rsite3[3], double r2[3], double v2[3] )
//------------------------------------------------------------------------------
/**
*
* @file anglesdr.cpp
*
* This function solves the problem of orbit determination using three optical sightings.
*
* @param az1 azimuth at t1 [rad]
* @param az2 azimuth at t2 [rad]
* @param az3 azimuth at t3 [rad]
* @param el1 elevation at t1 [rad]
* @param el2 elevation at t2 [rad]
* @param el3 elevation at t3 [rad]
* @param Mjd1 Modified julian date of t1
* @param Mjd2 Modified julian date of t2
* @param Mjd3 Modified julian date of t3
* @param rsite1 ijk site1 position vector [m]
* @param rsite2 ijk site2 position vector [m]
* @param rsite3 ijk site3 position vector [m]
* @param r ijk position vector at t2 [m]
* @param v ijk velocity vector at t2 [m/s]
*/ 
//------------------------------------------------------------------------------

extern double **eopdata;

void anglesdr ( double az1, double az2,double az3, double el1, double el2, double el3, double Mjd1, double Mjd2, double Mjd3, double rsite1[3], 
	double rsite2[3], double rsite3[3], double r2[3], double v2[3] ){

	cout << setprecision(16) << endl;

	double magr1in, magr2in, tol, pctchg, t1, t3, lon1, lon2, lon3, lat1, lat2, lat3, h1, h2, h3, M1[3][3], M2[3][3], M3[3][3], los1_m[3][1], los2_m[3][1],
	los3_m[3][1], los1_r[3][1], los2_r[3][1], los3_r[3][1], Mjd_UTC, UT1_UTC, TAI_UTC, x_pole, y_pole, UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC, Mjd_TT,
	Mjd_UT1, P[3][3], N[3][3], T[3][3], E[3][3], E1[3][3], E2[3][3], los1_r2[3][1], rsite1_m[3][1], rsite1_r[3][1], rsite2_m[3][1], rsite2_r[3][1],
	los2_r2[3][1], los3_r2[3][1], rsite3_m[3][1], rsite3_r[3][1], magr1old, magr2old, magrsite1, magrsite2, magrsite3, cc1, cc2, r3[3], magr1, magr2, a, 
	deltae32, f, g, magr1o, magr2o, rsite1_fin[3], rsite2_fin[3], rsite3_fin[3];
	char direct;
	complex<double> f1, f2, q1, q2, f1delr1, f2delr1, pf1pr1, pf2pr1, f1delr2,f2delr2,q3, pf1pr2, pf2pr2, delta, delta1, delta2, deltar1, deltar2;
	int ktr;
	
	magr1in = 1.1*R_Earth;
	magr2in = 1.11*R_Earth;
	direct  = 'y';
	
	tol    = 1e-8*R_Earth;
	pctchg = 0.005;
	
	t1 = (Mjd1 - Mjd2)*86400.0;
	t3 = (Mjd3 - Mjd2)*86400.0;
	
	double los1[3] = {cos(el1)*sin(az1), cos(el1)*cos(az1), sin(el1)};
	double los2[3] = {cos(el2)*sin(az2), cos(el2)*cos(az2), sin(el2)};
	double los3[3] = {cos(el3)*sin(az3), cos(el3)*cos(az3), sin(el3)};
	
	Geodetic(rsite1, &lon1, &lat1, &h1);
	Geodetic(rsite2, &lon2, &lat2, &h2);
	Geodetic(rsite3, &lon3, &lat3, &h3);
	
	LTC(lon1, lat1, M1);
	LTC(lon2, lat2, M2);
	LTC(lon3, lat3, M3);
	
	// body-fixed system
	transpuesta(M1, 3, 3);
	vectorAMatriz(los1, los1_m);
	vectorAMatriz(los2, los2_m);
	vectorAMatriz(los3, los3_m);
	
	mult3x1(M1, los1_m, los1_r);
	mult3x1(M1, los2_m, los2_r);
	mult3x1(M1, los3_m, los3_r);
	
	// mean of date system (J2000)
	Mjd_UTC = Mjd1;
	IERS(eopdata, Mjd_UTC, &UT1_UTC, &TAI_UTC, &x_pole, &y_pole);
	timediff(UT1_UTC, TAI_UTC, &UT1_TAI, &UTC_GPS, &UT1_GPS, &TT_UTC, &GPS_UTC);
	Mjd_TT = Mjd_UTC + TT_UTC/86400;
	Mjd_UT1 = Mjd_TT + (UT1_UTC-TT_UTC)/86400.0;
	
	PrecMatrix(MJD_J2000,Mjd_TT, P);
	NutMatrix(Mjd_TT, N);
	mult3x3(N, P, T);
	PoleMatrix(x_pole,y_pole, E1);
	GHAMatrix(Mjd_UT1, E2);
	mult3Matrices(E1, E2, T, E);
	
	transpuesta(E, 3, 3);
	mult3x1(E, los1_r, los1_r2);
	
	vectorAMatriz(rsite1, rsite1_m);
	mult3x1(E, rsite1_m, rsite1_r);
	
	Mjd_UTC = Mjd2;
	IERS(eopdata, Mjd_UTC, &UT1_UTC, &TAI_UTC, &x_pole, &y_pole);
	timediff(UT1_UTC, TAI_UTC, &UT1_TAI, &UTC_GPS, &UT1_GPS, &TT_UTC, &GPS_UTC);
	Mjd_TT = Mjd_UTC + TT_UTC/86400;
	Mjd_UT1 = Mjd_TT + (UT1_UTC-TT_UTC)/86400.0;
	
	PrecMatrix(MJD_J2000,Mjd_TT, P);
	NutMatrix(Mjd_TT, N);
	mult3x3(N, P, T);
	PoleMatrix(x_pole,y_pole, E1);
	GHAMatrix(Mjd_UT1, E2);
	mult3Matrices(E1, E2, T, E);
	
	transpuesta(E, 3, 3);
	mult3x1(E, los2_r, los2_r2);
	
	vectorAMatriz(rsite2, rsite2_m);
	mult3x1(E, rsite2_m, rsite2_r);
	
	Mjd_UTC = Mjd3;
	IERS(eopdata, Mjd_UTC, &UT1_UTC, &TAI_UTC, &x_pole, &y_pole);
	timediff(UT1_UTC, TAI_UTC, &UT1_TAI, &UTC_GPS, &UT1_GPS, &TT_UTC, &GPS_UTC);
	Mjd_TT = Mjd_UTC + TT_UTC/86400;
	Mjd_UT1 = Mjd_TT + (UT1_UTC-TT_UTC)/86400.0;
	
	PrecMatrix(MJD_J2000,Mjd_TT, P);
	NutMatrix(Mjd_TT, N);
	mult3x3(N, P, T);
	PoleMatrix(x_pole,y_pole, E1);
	GHAMatrix(Mjd_UT1, E2);
	mult3Matrices(E1, E2, T, E);
	
	transpuesta(E, 3, 3);
	mult3x1(E, los3_r, los3_r2);
	
	vectorAMatriz(rsite3, rsite3_m);
	mult3x1(E, rsite3_m, rsite3_r);
	
	matrizAVector(rsite1_r, rsite1_fin);
	matrizAVector(rsite2_r, rsite2_fin);
	matrizAVector(rsite3_r, rsite3_fin);
	
	
	matrizAVector(los1_r2, los1);
	matrizAVector(los2_r2, los2);
	matrizAVector(los3_r2, los3);
	
	magr1old  = 99999999.9;
	magr2old  = 99999999.9;
	magrsite1 = norm(rsite1_fin, 3);
	magrsite2 = norm(rsite2_fin, 3);
	magrsite3 = norm(rsite3_fin, 3);
	
	
	cc1 = 2.0*dot(los1, 3, rsite1_fin, 3);
	cc2 = 2.0*dot(los2, 3, rsite2_fin, 3);
	ktr = 0;
	
		//Comprobado hasta aqui
	
	while (fabs(magr1in-magr1old) > tol | fabs(magr2in-magr2old) > tol){
		ktr = ktr + 1;
		doubler( cc1,cc2,magrsite1,magrsite2,magr1in,magr2in,los1,los2,los3,rsite1_fin,rsite2_fin,rsite3_fin,t1,t3,direct, r2,r3,&f1,&f2,&q1,&magr1,&magr2,&a,&deltae32);
		
		f  = 1.0 - a/magr2*(1.0-cos(deltae32));
		g  = t3 - sqrt(a*a*a/GM_Earth)*(deltae32-sin(deltae32));
		for(int i=0; i<3; i++){
			v2[i] = (r3[i] -f*r2[i])/g;
		}
		
		magr1o = magr1in;
		magr1in = (1.0+pctchg)*magr1in;
		deltar1 = pctchg*magr1in;
	    doubler( cc1,cc2,magrsite1,magrsite2,magr1in,magr2in,los1,los2,los3,rsite1_fin,rsite2_fin,rsite3_fin,t1,t3,direct, r2,r3,&f1delr1,&f2delr1,&q2,&magr1,&magr2,&a,&deltae32);
		
		pf1pr1 = (f1delr1-f1)/deltar1;
		pf2pr1 = (f2delr1-f2)/deltar1;
		
		magr1in = magr1o;
		deltar1 = pctchg*magr1in;
		magr2o = magr2in;
		magr2in = (1.0+pctchg)*magr2in;
		deltar2 = pctchg*magr2in;
		doubler( cc1,cc2,magrsite1,magrsite2,magr1in,magr2in, los1,los2,los3,rsite1_fin,rsite2_fin,rsite3_fin,t1,t3,direct, r2,r3,&f1delr2,&f2delr2,&q3,&magr1,&magr2,&a,&deltae32);
		
		pf1pr2 = (f1delr2-f1)/deltar2;
		pf2pr2 = (f2delr2-f2)/deltar2;
		
		magr2in = magr2o;
		deltar2 = pctchg*magr2in;
		
		delta  = pf1pr1*pf2pr2 - pf2pr1*pf1pr2;
		delta1 = pf2pr2*f1 - pf1pr2*f2;
		delta2 = pf1pr1*f2 - pf2pr1*f1;
		
		deltar1 = -delta1/delta;
		deltar2 = -delta2/delta;
		
		magr1old = magr1in;
		magr2old = magr2in;
		
		magr1in = magr1in + fabs(deltar1);
		magr2in = magr2in + fabs(deltar2);
		
	}
	
	doubler( cc1,cc2,magrsite1,magrsite2,magr1in,magr2in,los1,los2,los3,rsite1_fin,rsite2_fin,rsite3_fin,t1,t3,direct, r2,r3,&f1,&f2,&q1,&magr1,&magr2,&a,&deltae32);
	
	f  = 1.0 - a/magr2*(1.0-cos(deltae32));
	g  = t3 - sqrt(a*a*a/GM_Earth)*(deltae32-sin(deltae32));
	for(int i=0; i<3; i++){
		v2[i] = (r3[i] -f*r2[i])/g;
	}
	
}