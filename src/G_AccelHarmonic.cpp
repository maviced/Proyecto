#include <cmath>
#include <iostream>
#include "../include/G_AccelHarmonic.h"
#include "../include/AccelHarmonic.h"
#include "../include/matriz.h"

//------------------------------------------------------------------------------
// void G_AccelHarmonic(double r[3][1], double U[3][3], int n_max, int m_max,double G[3][3])
//------------------------------------------------------------------------------
/**
*
* @file G_AccelHarmonic.cpp
*
* Computes the gradient of the Earth's harmonic gravity field 
*
* @param Mjd_UT Modified Julian Date (Universal Time)
* @param r Satellite position vector in the true-of-date system
* @param n_max Gravity model degree
* @param m_max Gravity model order
* @param G Gradient (G=da/dr) in the true-of-date system
*/ 
//------------------------------------------------------------------------------

using namespace std;

void G_AccelHarmonic(double r[3][1], double U[3][3], int n_max, int m_max,double G[3][3]){
	
	double d, dr[3][1], da[3][1], da1[3][1], da2[3][1], aux1[3][1], aux2[3][1];
	
	d = 1.0;
	
	for(int i=0; i<3; i++){
		// Set offset in i-th component of the position vector
		dr[0][0] = dr[1][0] = dr[2][0] = 0.0;
		dr[i][0] = d;
		
		// Acceleration difference
		for (int j = 0; j < 3; j++) {
			aux1[j][0] = r[j][0] + dr[j][0]/2;
			aux2[j][0] = r[j][0] - dr[j][0]/2;
		}
		AccelHarmonic( aux1,U, n_max, m_max, da1 );
		transpuesta(U, 3, 3);
		AccelHarmonic( aux2,U, n_max, m_max, da2 );
		for (int j = 0; j < 3; j++) {
			da[j][0] = da1[j][0] - da2[j][0];
		}
		// Derivative with respect to i-th axis
		//G(:,i) = da/d;   
		for (int j = 0; j < 3; j++) {
			G[j][i] = da[j][0]/d;
		}
	}
}