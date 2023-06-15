#include <cmath>
#include "../include/NutAngles.h"
#include "../include/EqnEquinox.h"
#include "../include/MeanObliquity.h"

//------------------------------------------------------------------------------
// double EqnEquinox(double Mjd_TT)
//------------------------------------------------------------------------------
/**
*
* @file EqnEquinox.cpp
*
* Computation of the equation of the equinoxes
*
* @param Mjd_TT Modified Julian Date (Terrestrial Time)
* @return Equation of the equinoxes
*
* @note The equation of the equinoxes dpsi*cos(eps) is the right ascension of
*   the mean equinox referred to the true equator and equinox and is equal
*   to the difference between apparent and mean sidereal time.
*/ 
//------------------------------------------------------------------------------

double EqnEquinox(double Mjd_TT){
	
	double dpsi, deps, EqE;
	
	// Nutation in longitude and obliquity
	NutAngles (Mjd_TT, &dpsi, &deps);

	// Equation of the equinoxes
	EqE = dpsi * cos ( MeanObliquity(Mjd_TT) );
	
	return EqE;
}