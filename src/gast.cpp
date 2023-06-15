#include <cmath>
#include "../include/gast.h"
#include "../include/gmst.h"
#include "../include/EqnEquinox.h"
#include "../include/SAT_Const.h"

//------------------------------------------------------------------------------
// double gast(double Mjd_UT1)
//------------------------------------------------------------------------------
/**
*
* @file gast.cpp
*
* Greenwich Apparent Sidereal Time
*
* @param Mjd_UT1 Modified Julian Date UT1
* @return GAST in [rad]
*/ 
//------------------------------------------------------------------------------

double gast(double Mjd_UT1){
	
	return fmod ( gmst(Mjd_UT1) + EqnEquinox(Mjd_UT1), 2*pi );
}