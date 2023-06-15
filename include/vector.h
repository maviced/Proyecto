//$Header$
//------------------------------------------------------------------------------
// vector
//------------------------------------------------------------------------------
// EKF_GEOS3: Initial Orbit Determination using Double-R-Iteration and Extended 
// Kalman Filter methods
//
//
// Last modified:   2015/08/12   M. Mahooti
//
/**
*
* @file vector.h
*
* Provides basic operations to manipulate static vectors of different sizes.
*
*/
//------------------------------------------------------------------------------

#ifndef _NORM_
#define _NORM_

double norm(double v[], int n);

#endif

#ifndef _DOT_
#define _DOT_

double dot(double v[], int n, double w[], int m);

#endif

#ifndef _CROSS_
#define _CROSS_

void cross(double v[], int n, double w[], int m, double s[], int &i);

#endif

#ifndef _vectoresIguales_
#define _vectoresIguales_

bool vectoresIguales(double v1[3], double v2[3], int n);

#endif

#ifndef _vIgualesMenosPrecision_
#define _vIgualesMenosPrecision_

bool vIgualesMenosPrecision(double v1[3], double v2[3], int n);

#endif

