#include <cmath>
#include <iostream>
#include "../include/VarEqn.h"
#include "../include/AccelHarmonic.h"
#include "../include/G_AccelHarmonic.h"
#include "../include/SAT_Const.h"
#include "../include/IERS.h"
#include "../include/timediff.h"
#include "../include/globales.h"
#include "../include/PrecMatrix.h"
#include "../include/NutMatrix.h"
#include "../include/matriz.h"
#include "../include/PoleMatrix.h"
#include "../include/GHAMatrix.h"
#include "../include/matriz.h"

//------------------------------------------------------------------------------
// void VarEqn(double x, double yPhi[42][1], double yPhip[42][1])
//------------------------------------------------------------------------------
/**
*
* @file VarEqn.cpp
*
* Computes the variational equations, i.e. the derivative of the state vector
*   and the state transition matrix
*
* @param x Time since epoch in [s]
* @param yPhi (6+36)-dim vector comprising the state vector (y) and the state transition matrix (Phi) in column wise storage order
* @param yPhip Derivative of yPhi
*/ 
//------------------------------------------------------------------------------

using namespace std;

void VarEqn(double x, double yPhi[42][1], double yPhip[42][1]){
	
	double UT1_UTC, TAI_UTC, x_pole, y_pole, UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC, Mjd_UT1, P[3][3], N[3][3], T[3][3], E1[3][3], E2[3][3], E[3][3],
	Phi[6][6], a[3][1], G[3][3], dfdy[6][6], Phip[6][6];
	
	a[0][0]=0.0;
	a[1][0]=0.0;
	a[2][0]=0.0;
	
	IERS(eopdata, AuxParam.Mjd_TT + x/86400, &UT1_UTC, &TAI_UTC, &x_pole, &y_pole);
	timediff(UT1_UTC, TAI_UTC, &UT1_TAI, &UTC_GPS, &UT1_GPS, &TT_UTC, &GPS_UTC);
	Mjd_UT1 = AuxParam.Mjd_TT + (UT1_UTC-TT_UTC)/86400.0;
	
	// Transformation matrix
	PrecMatrix(MJD_J2000,AuxParam.Mjd_TT + x/86400, P);
	NutMatrix(AuxParam.Mjd_TT + x/86400, N);
	mult3x3(N, P, T);
	PoleMatrix(x_pole,y_pole, E1);
	GHAMatrix(Mjd_UT1, E2);
	mult3Matrices(E1, E2, T, E);
	
	
	// State vector components
	double r[3][1] = {{yPhi[0][0]}, {yPhi[1][0]}, {yPhi[2][0]}};
	double v[3][1] = {{yPhi[3][0]}, {yPhi[4][0]}, {yPhi[5][0]}};
	
	AccelHarmonic ( r, E, AuxParam.n, AuxParam.m , a);
	transpuesta(E, 3, 3);
	G_AccelHarmonic( r, E, AuxParam.n, AuxParam.m , G);
	
	// State transition matrix
	for (int j=1; j<7; j++){
		for(int i=0; i<6; i++){
			Phi[i][j-1] = yPhi[6*j+i][0];
		}
	}
	
	
	// Acceleration and gradient
	//AccelHarmonic ( r, E, AuxParam.n, AuxParam.m , a);
	//G_AccelHarmonic ( r, E, AuxParam.n, AuxParam.m , G);

	for (int i=0; i<3; i++){
		for (int j=0; j<3; j++){
			dfdy[i][j] = 0.0;                 // dv/dr(i,j)
			dfdy[i+3][j] = G[i][j];            // da/dr(i,j)
			if ( i==j ){
				dfdy[i][j+3] = 1.0;
			}else{
				dfdy[i][j+3] = 0.0;             // dv/dv(i,j)
			}
			dfdy[i+3][j+3] = 0.0;             // da/dv(i,j)
		}
	}
	
	mult6x6(dfdy, Phi, Phip);
	
	// Derivative of combined state vector and state transition matrix
	for (int i=0; i<3; i++){
		yPhip[i][0]   = v[i][0];                 // dr/dt(i)
		yPhip[i+3][0] = a[i][0];                 // dv/dt(i)
	}
	
	for (int i=0; i<6; i++){
		for (int j=1; j<7; j++){
			yPhip[6*j+i][0] = Phip[i][j-1];     // dPhi/dt(i,j)
		}
	}
	
}