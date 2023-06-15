#include <cmath>
#include "../include/Mjday.h"

//------------------------------------------------------------------------------
// double Mjday (int yr, int mon, int day, int hr, int min, double sec)
//------------------------------------------------------------------------------
/**
*
* @file Mjday.cpp
*
* Modified Julian Day (MJD) calculation based on date and time components.
*
* @param yr year
* @param mon month
* @param day day
* @param hr universal time hour
* @param min universal time min
* @param sec universal time sec
* @return Modified julian date
*/ 
//------------------------------------------------------------------------------

double Mjday (int yr, int mon, int day, int hr, int min, double sec){
	
	double jd = 367.0 * yr - floor( (7 * (yr + floor( (mon + 9) / 12.0) ) ) * 0.25 ) + floor( 275 * mon / 9.0 ) + day + 1721013.5 + ( (sec/60.0 + min ) / 60.0 + hr ) / 24.0;
	
	double Mjd = jd-2400000.5;
	
	return Mjd;
}


