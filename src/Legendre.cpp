#include <cmath>
#include <iostream>
#include <iomanip>
#include "../include/TimeUpdate.h"

//------------------------------------------------------------------------------
// void Legendre(int n,int m,double fi, double pnm[362][362], double dpnm[362][362])
//------------------------------------------------------------------------------
/**
*
* @file Legendre.cpp
*
* Implements the computation of Legendre polynomials and their derivatives using a recurrence relation.
*
* @param n The maximum degree of the Legendre polynomials.
* @param m The maximum order of the Legendre polynomials.
* @param fi A double representing an angle [rad].
* @param pnm An array of size [362][362] to store the Legendre polynomials.
* @param dpnm An array of size [362][362] to store the derivatives of the Legendre polynomials.
*/ 
//------------------------------------------------------------------------------

using namespace std;

void Legendre(int n,int m,double fi, double pnm[362][362], double dpnm[362][362]){
	
	int j, k;
	
	pnm[0][0]=1;
	dpnm[0][0]=0;
	pnm[1][1]=sqrt(3)*cos(fi);
	dpnm[1][1]=-sqrt(3)*sin(fi);
	
	// diagonal coefficients
	for (int i=2; i<n+1; i++){  
		pnm[i][i]= sqrt((2.0*i+1)/(2.0*i))*cos(fi)*pnm[i-1][i-1];
	}
	for (int i=2; i<n+1; i++){
		dpnm[i][i]= sqrt((2.0*i+1)/(2.0*i))*((cos(fi)*dpnm[i-1][i-1])-(sin(fi)*pnm[i-1][i-1]));
	}
	
	// horizontal first step coefficients
	for (int i=1; i<n+1; i++){ 
		pnm[i][i-1]= sqrt(2.0*i+1)*sin(fi)*pnm[i-1][i-1];
	}
	for (int i=1; i<n+1; i++){ 
		dpnm[i][i-1]= sqrt(2.0*i+1)*((cos(fi)*pnm[i-1][i-1])+(sin(fi)*dpnm[i-1][i-1]));
	}
	
	// horizontal second step coefficients
	j=0;
	k=2;
	while(1){
		for (int i=k; i<n+1; i++){     
			pnm[i][j]=sqrt((2.0*i+1)/((i-j)*(i+j)))*((sqrt(2.0*i-1)*sin(fi)*pnm[i-1][j])-(sqrt(((i+j-1.0)*(i-j-1))/(2.0*i-3))*pnm[i-2][j]));
		}
		j = j+1;
		k = k+1;
		if (j>m){
			break;
		}
	}
	j = 0;
	k = 2;
	while(1){
		for (int i=k; i<n+1; i++){       
			dpnm[i][j]=sqrt((2.0*i+1)/((i-j)*(i+j)))*((sqrt(2.0*i-1)*sin(fi)*dpnm[i-1][j])+(sqrt(2.0*i-1)*cos(fi)*pnm[i-1][j])-(sqrt(((i+j-1.0)*(i-j-1))/(2.0*i-3))*dpnm[i-2][j]));
		}
		j = j+1;
		k = k+1;
		if (j>m){
			break;
		}
	}
		
}