#include <cmath>
#include "../include/timediff.h"

//------------------------------------------------------------------------------
// void timediff(double UT1_UTC, double TAI_UTC, double *UT1_TAI, double *UTC_GPS, double *UT1_GPS, double *TT_UTC, double *GPS_UTC)
//------------------------------------------------------------------------------
/**
*
* @file timediff.cpp
*
* Time differences [s]
*
* @param UT1_UTC the difference between UT1 and UTC [s]
* @param TAI_UTC the difference between TAI and UTC [s]
* @param UT1_TAI pointer to store the difference between UT1 and TAI [s]
* @param UTC_GPS pointer to store the difference between UTC and GPS [s]
* @param UT1_GPS pointer to store the difference between UT1 and GPS [s]
* @param TT_UTC pointer to store the difference between TT and UTC [s]
* @param GPS_UTC pointer to store the difference between GPS and UTC [s]
*/ 
//------------------------------------------------------------------------------

void timediff(double UT1_UTC, double TAI_UTC, double *UT1_TAI, double *UTC_GPS, double *UT1_GPS, double *TT_UTC, double *GPS_UTC){

	double TT_TAI, GPS_TAI, TT_GPS, TAI_GPS, UTC_TAI;

	TT_TAI  = +32.184;          // TT-TAI time difference [s]

	GPS_TAI = -19.0;            // GPS-TAI time difference [s]

	TT_GPS  =  TT_TAI-GPS_TAI;  // TT-GPS time difference [s]

	TAI_GPS = -GPS_TAI;         // TAI-GPS time difference [s]

	*UT1_TAI = UT1_UTC-TAI_UTC;  // UT1-TAI time difference [s]

	UTC_TAI = -TAI_UTC;         // UTC-TAI time difference [s]
	  
	*UTC_GPS = UTC_TAI-GPS_TAI;  // UTC_GPS time difference [s]

	*UT1_GPS = *UT1_TAI-GPS_TAI;  // UT1-GPS time difference [s]

	*TT_UTC  = TT_TAI-UTC_TAI;   //  TT-UTC time difference [s]

	*GPS_UTC = GPS_TAI-UTC_TAI;  // GPS-UTC time difference [s]


}