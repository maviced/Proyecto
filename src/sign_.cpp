#include <cmath>
#include "../include/sign_.h"

//------------------------------------------------------------------------------
// double sign_(double a,double b)
//------------------------------------------------------------------------------
/**
*
* @file sign_.cpp
*
* Returns absolute value of a with sign of b.
*
* @param a The value whose signed magnitude needs to be calculated.
* @param b Determines the sign of the result.
* @return Absolute value of a with sign of b.
*/ 
//------------------------------------------------------------------------------

double sign_(double a,double b){

	double result;

	if (b>=0.0){
		result = fabs(a);
	}else{
		result = - fabs(a);
	}
	return result;
}