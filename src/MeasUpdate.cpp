#include <cmath>
#include "../include/TimeUpdate.h"
#include "../include/matriz.h"

//------------------------------------------------------------------------------
// void MeasUpdate(double x[6][1], double  z, double g, double s,double G[1][6], double P[6][6],int n, double K[6][1], double x2[6][1],double  P2[6][6])
//------------------------------------------------------------------------------
/**
*
* @file MeasUpdate.cpp
*
* Implements a measurement update step for a Kalman filter
*
* @param x The state vector.
* @param z The measured value.
* @param g The expected measurement value based on the state vector.
* @param s The standard deviation of the measurement noise.
* @param G The measurement sensitivity matrix.
* @param P The covariance matrix.
* @param n The dimensionality of the state vector (in this case, 6).
* @param K The Kalman gain.
* @param x2 The updated state vector.
* @param P2 The updated covariance matrix.
*/ 
//------------------------------------------------------------------------------

void MeasUpdate(double x[6][1], double  z, double g, double s,double G[1][6], double P[6][6],int n, double K[6][1], double x2[6][1],double  P2[6][6]){
	double Inv_W, result1[1][6], G_trans[6][1], result2, aux[6][1], aux2[6][6], result3[6][6], result4[6][6];
	
	Inv_W = s*s;
	
	// Kalman gain
	//K = P*G'*inv(Inv_W+G*P*G');
	multMatrices1x6_6x6(G, P, result1);
	transpuesta1x6(G, G_trans);
	result2 = 1/(multMatrices1x6_6x1(result1, G_trans) + Inv_W);
	multMatrices6x6_6x1(P, G_trans, K);
	for (int i =0; i<6; i++){
		K[i][0] = K[i][0] *result2;
	}
	
	// State update
	//x = x + K*(z-g);
	for (int i =0; i<6; i++){
		aux[i][0] = K[i][0] *(z-g);
	}
	for(int j=0; j<6; j++){
		x2[j][0] = x[j][0] +aux[j][0];
	}
	
	// Covariance update
	//P = (eye(n)-K*G)*P;
	eye6x6(aux2);
	multMatrices6x1_1x6(K, G, result3);
	for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            result4[i][j] = aux2[i][j] - result3[i][j];
        }
    }
	multMatrices6x6(result4, P, P2);
	
	
}