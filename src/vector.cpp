#include <cmath>
#include "../include/vector.h"

/**
*
* @file vector.cpp
*
* Provides basic operations to manipulate static vectors of different sizes.
*
*/


//------------------------------------------------------------------------------
// double norm(double v[], int n)
//------------------------------------------------------------------------------
/**
*
* Calculates the norm (magnitude) of a vector.
*
* @param v The input vector.
* @param n The number of elements in the vector.
* @return The norm of the vector.
* @throws An exception "Empty vector" if the vector is empty (n <= 0).
*/ 
//------------------------------------------------------------------------------

double norm(double v[], int n)
{
	double suma=0;
	int i;
	
	if(n <= 0)
		throw "Empty vector";
	
	for(i=0; i<n; i++)
		suma += v[i]*v[i];
	
	return(sqrt(suma));
}


//------------------------------------------------------------------------------
// double dot(double v[], int n, double w[], int m)
//------------------------------------------------------------------------------
/**
*
* Calculates the dot product of two vectors.
*
* @param v The first input vector.
* @param n The number of elements in the first vector.
* @param w The second input vector.
* @param m The number of elements in the second vector.
* @return The dot product of the two vectors.
* @throws An exception "Empty vector" if either of the vectors is empty (n <= 0 or m <= 0).
* @throws An exception "Different dimensions" if the dimensions of the two vectors are different (n != m).
*/ 
//------------------------------------------------------------------------------

double dot(double v[], int n, double w[], int m)
{
	double suma=0;
	int i;
	
	if( (n <= 0) || (m <=0))
		throw "Empty vector";
	
	if(n != m)
		throw "Different dimensions";
	
	for(i=0; i<n; i++)
		suma += v[i]*w[i];
	
	return(suma);
}


//------------------------------------------------------------------------------
// void cross(double v[], int n, double w[], int m, double s[], int &i)
//------------------------------------------------------------------------------
/**
*
* Calculates the cross product of two vectors.
*
* @param v The first input vector.
* @param n The number of elements in the first vector.
* @param w The second input vector.
* @param m The number of elements in the second vector.
* @param s The resulting cross product vector.
* @param i The number of elements in the resulting cross product vector.
* @throws An exception "No valid dimensions" if the dimensions of the input vectors are not valid (n != 3 or m != 3).
*/ 
//------------------------------------------------------------------------------


void cross(double v[], int n, double w[], int m, double s[], int &i)
{
	
	if ( (n != 3) || (m != 3))
		throw "No valid dimensions";
	
	i=n;
	s[0] = v[1]*w[2]-v[2]*w[1];
	s[1] = v[2]*w[0]-v[0]*w[2];
	s[2] = v[0]*w[1]-v[1]*w[0];
	
	
}

//------------------------------------------------------------------------------
// bool vectoresIguales(double v1[3], double v2[3], int n)
//------------------------------------------------------------------------------
/**
*
* Checks if two vectors are equal within a specified tolerance.
*
* @param v1 The first input vector.
* @param v2 The second input vector.
* @param n The number of elements in the vectors.
* @return True if the vectors are equal within the specified tolerance, false otherwise.
*/ 
//------------------------------------------------------------------------------

bool vectoresIguales(double v1[3], double v2[3], int n) {
    
    for (int i = 0; i < n; i++) {
	
		if (fabs(v1[i] - v2[i]) > pow(10,-9)) {
			return false;
		}
	
    }
    return true;
}

//------------------------------------------------------------------------------
// bool vIgualesMenosPrecision(double v1[3], double v2[3], int n)
//------------------------------------------------------------------------------
/**
*
* Checks if two vectors are equal within a slightly lower precision tolerance.
*
* @param v1 The first input vector.
* @param v2 The second input vector.
* @param n The number of elements in the vectors.
* @return True if the vectors are equal within the slightly lower precision tolerance, false otherwise.
*/ 
//------------------------------------------------------------------------------

bool vIgualesMenosPrecision(double v1[3], double v2[3], int n){
	for (int i = 0; i < n; i++) {
	
		if (fabs(v1[i] - v2[i]) > pow(10,-6)) {
			return false;
		}
	
    }
    return true;
}
