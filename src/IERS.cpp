#include <cmath>
#include <iostream>
#include "../include/IERS.h"
#include "../include/SAT_Const.h"

//------------------------------------------------------------------------------
// void IERS(double **eop, double Mjd_TT, double *UT1_UTC, double *TAI_UTC, double *x_pole, double *y_pole)
//------------------------------------------------------------------------------
/**
*
* @file IERS.cpp
*
* Management of IERS time and polar motion data
*
* @param eop Earth Orientation Parameters (EOP) data
* @param Mjd_TT Modified Julian Date (MJD) for Terrestrial Time (TT)
* @param UT1_UTC The UT1-UTC time difference in seconds
* @param TAI_UTC The TAI-UTC time difference in seconds
* @param x_pole The pole coordinate in radians (converted from arcseconds)
* @param y_pole The pole coordinate in radians (converted from arcseconds)
*/ 
//------------------------------------------------------------------------------


using namespace std;

void IERS(double **eop, double Mjd_TT, double *UT1_UTC, double *TAI_UTC, double *x_pole, double *y_pole){
	
	double Arcs, v[13];
	int mj, nop;
	
	Arcs = 3600.0*180.0/pi;  // Arcseconds per radian
	
	mj = (round(Mjd_TT));
	nop = 19716;
	
	/*for(int i=0; i<nop; i++){
		if (mj==eop[i][3]){
			for(int j=0;j<13;j++){
				v[j] = eop[i][j];
			}
			break;
		}
	}*/
	
	for(int i=0; i<nop; i++){
		if (mj==eop[3][i]){
			for(int j=0;j<13;j++){
				v[j] = eop[j][i];
			}
			break;
		}
	}
	
	
	
	// Setting of IERS Earth rotation parameters
	// (UT1-UTC [s], TAI-UTC [s], x ["], y ["])
	*UT1_UTC = v[6];     // UT1-UTC time difference [s]
	*TAI_UTC = v[12];     // TAI-UTC time difference [s]
	*x_pole  = v[4]/Arcs; // Pole coordinate [rad]
	*y_pole  = v[5]/Arcs; // Pole coordinate [rad]
}
	