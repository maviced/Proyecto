#include <cmath>
#include "../include/Frac.h"

//------------------------------------------------------------------------------
// double Frac (double x)
//------------------------------------------------------------------------------
/**
*
* @file seebatt.cpp
*
* Fractional part of a number (y=x-[x])
*
* @param x The input number for which the fractional part needs to be calculated.
* @return The fractional part of the input number.
*/ 
//------------------------------------------------------------------------------


double Frac (double x){
	
	double res = x - floor (x);
	
	return res;
}
