#ifndef _globales_
#define _globales_

extern double **eopdata;
extern double pnm[362][362], dpnm[362][362], Cnm[362][362], Snm[362][362];
typedef struct{
	double Mjd_TT;
    int n;
	int m;
    double Mjd_UTC;
} Param;
extern Param AuxParam;

#endif

// = 4.974611324287046e+04
// = 4.974611253472231e+04