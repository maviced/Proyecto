#include <cmath>
#include "../include/TimeUpdate.h"
#include "../include/matriz.h"

//------------------------------------------------------------------------------
// void TimeUpdate(double P[6][6], double Phi[6][6], double Qdt)
//------------------------------------------------------------------------------
/**
*
* @file TimeUpdate.cpp
*
* Implements a time update step for a Kalman filter.
*
* @param P Covariance matrix.
* @param Phi State transition matrix.
* @param Qdt The process noise covariance multiplied by the time step.
*/ 
//------------------------------------------------------------------------------

void TimeUpdate(double P[6][6], double Phi[6][6], double Qdt){
	
	double Phi_aux[6][6], P_aux[6][6];
	
	copiarMatriz6x6(Phi, Phi_aux);
	transpuesta6x6(Phi,6,6);
	mult3Matrices6x6(Phi_aux, P, Phi, P_aux);
	
	copiarMatriz6x6(P_aux, P);
	
	for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            P[i][j] = P[i][j] + Qdt;
        }
    }
	
}