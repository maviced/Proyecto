#include <cmath>
#include <complex>
#include <iostream>
#include <iomanip>
#include "../include/doubler.h"
#include "../include/SAT_Const.h"
#include "../include/vector.h"

using namespace std;

//------------------------------------------------------------------------------
// void doubler(double cc1,double cc2,double magrsite1,double magrsite2,double magr1in,double magr2in,double los1[3],
// double los2[3],double los3[3],double rsite1[3],double rsite2[3],double rsite3[3],double t1,double t3,char direct,
// double r2[3],double r3[3],complex<double> *f1,complex<double> *f2,complex<double> *q1,double *magr1,double *magr2,double *a,double *deltae32)
//------------------------------------------------------------------------------
/**
*
* @file doubler.cpp
*
* Performs various calculations related to double-r observations.
*
* @param cc1 constant value
* @param cc2 constant value
* @param magrsite1 magnitude of rsite1 vector
* @param magrsite2 magnitude of rsite2 vector
* @param magr1in magnitude of magr1in vector
* @param magr2in magnitude of magr2in vector
* @param los1 array of size 3 containing los1 values
* @param los2 array of size 3 containing los2 values
* @param los3 array of size 3 containing los3 values
* @param rsite1 array of size 3 containing rsite1 values
* @param rsite2 array of size 3 containing rsite2 values
* @param rsite3 array of size 3 containing rsite3 values
* @param t1 value of t1
* @param t3 value of t3
* @param direct character indicating the direction
* @param r2 array of size 3 to store r2 values
* @param r3 array of size 3 to store r3 values
* @param f1 pointer to store f1 value
* @param f2 pointer to store f2 value
* @param q1 pointer to store q1 value
* @param magr1 pointer to store magr1 value
* @param magr2 pointer to store magr2 value
* @param a pointer to store a value
* @param deltae32 pointer to store deltae32 value
*/ 
//------------------------------------------------------------------------------

void doubler(double cc1,double cc2,double magrsite1,double magrsite2,double magr1in,double magr2in,double los1[3],
 double los2[3],double los3[3],double rsite1[3],double rsite2[3],double rsite3[3],double t1,double t3,char direct,
 double r2[3],double r3[3],complex<double> *f1,complex<double> *f2,complex<double> *q1,double *magr1,double *magr2,double *a,double *deltae32){
	 
	double rho1, rho2, r1[3], aux[3], w[3], rho3, magr3, cosdv21, sindv21, dv21, cosdv31, sindv31, dv31, cosdv32, sindv32, dv32, c1, c3, p,
	ecosv1, ecosv2, ecosv3, esinv2, e, s, c, sinde32, cosde32, sinde21, cosde21, deltam32, deltam12, deltae21, sindh32, sindh21, deltah32,
	deltah21;
	complex<double> n;
	int i;
	
	rho1 = (-cc1 + sqrt(cc1*cc1-4*(magrsite1*magrsite1-magr1in*magr1in))) / 2.0;
	rho2 = (-cc2 + sqrt(cc2*cc2-4*(magrsite2*magrsite2-magr2in*magr2in))) / 2.0;
	
	
	for(int i=0; i<3; i++){
		r1[i] = rho1*los1[i] + rsite1[i];
		r2[i] = rho2*los2[i] + rsite2[i];
	}
	
	*magr1 = norm(r1, 3);
	*magr2 = norm(r2, 3);
	
	
	if (direct == 'y'){
		cross(r1,3,r2,3,aux,i);
		for(int i=0; i<3; i++){
			w[i] = aux[i] /(*magr1 * *magr2);
		}
	}else{
		cross(r1,3,r2,3,aux,i);
		for(int i=0; i<3; i++){
			w[i] = -aux[i] /(*magr1 * *magr2);
		}
	}
	
	rho3 =  -dot(rsite3,3,w,3)/dot(los3,3,w,3);
	for(int i=0;i<3;i++){
		r3[i] = rho3*los3[i] + rsite3[i];
	}
	magr3 = norm(r3, 3);
	
	cosdv21 = dot(r2,3,r1,3)/(*magr2 * *magr1);
	cross(r2,3,r1,3,aux,i);
	sindv21 = norm(aux,3)/(*magr2 * *magr1);
	dv21 = atan2(sindv21,cosdv21);
	
	cosdv31 = dot(r3,3,r1,3)/(magr3* *magr1);
	sindv31 = sqrt(1.0 - cosdv31*cosdv31);
	dv31 = atan2(sindv31,cosdv31);
	
	cosdv32 = dot(r3,3,r2,3)/(magr3* *magr2);
	cross(r3,3,r2,3,aux,i);
	sindv32 = norm(aux,3)/(magr3* *magr2);
	dv32 = atan2(sindv32,cosdv32); 
	
	if (dv31 > pi){
		c1 = (*magr2*sindv32)/(*magr1*sindv31);
		c3 = (*magr2*sindv21)/(magr3*sindv31);
		p = (c1* *magr1+c3*magr3-*magr2)/(c1+c3-1);
	}else{
		c1 = (*magr1*sindv31)/(*magr2*sindv32);
		c3 = (*magr1*sindv21)/(magr3*sindv32);
		p = (c3*magr3-c1* *magr2+*magr1)/(-c1+c3+1);
	}
	
	ecosv1 = p/ *magr1-1;
	ecosv2 = p/ *magr2-1;
	ecosv3 = p/magr3-1;
	
	
	if (dv21!=pi){
		esinv2 = (-cosdv21*ecosv2+ecosv1)/sindv21;
	}else{
		esinv2 = (cosdv32*ecosv2-ecosv3)/sindv31;
	}
	
	e = sqrt(ecosv2*ecosv2+esinv2*esinv2);
	*a = p/(1-e*e);
	
	if (e*e < 0.99){
		n = sqrt(GM_Earth/(*a* *a* *a));
		
		s = *magr2/p*sqrt(1-e*e)*esinv2;
		c = *magr2/p*(e*e+ecosv2);
		
		sinde32 = magr3/sqrt(*a*p)*sindv32-magr3/p*(1-cosdv32)*s;
		cosde32 = 1-*magr2*magr3/(*a*p)*(1-cosdv32);
		*deltae32 = atan2(sinde32,cosde32);
		
		sinde21 = *magr1/sqrt(*a*p)*sindv21+*magr1/p*(1-cosdv21)*s;
		cosde21 = 1-*magr2* *magr1/(*a*p)*(1-cosdv21);
		deltae21 = atan2(sinde21,cosde21);
		
		deltam32 = *deltae32+2*s*(sin(*deltae32/2))*(sin(*deltae32/2))-c*sin(*deltae32);
		deltam12 = -deltae21+2*s*(sin(deltae21/2))*(sin(deltae21/2))+c*sin(deltae21);
	}else{
		n = sqrt(complex<double>(GM_Earth/-(*a* *a* *a)));
		
		s = *magr2/p*sqrt(e*e-1)*esinv2;
		c = *magr2/p*(e*e+ecosv2);
		
		sindh32 = magr3/sqrt(-*a*p)*sindv32-magr3/p*(1-cosdv32)*s;
		sindh21 = *magr1/sqrt(-*a*p)*sindv21+*magr1/p*(1-cosdv21)*s;
		
		deltah32 = log( sindh32 + sqrt(sindh32*sindh32 +1) );
		deltah21 = log( sindh21 + sqrt(sindh21*sindh21 +1) );
		
		deltam32 = -deltah32+2*s*(sinh(deltah32/2))*(sinh(deltah32/2))+c*sinh(deltah32);
		deltam12 = deltah21+2*s*(sinh(deltah21/2))*(sinh(deltah21/2))-c*sinh(deltah21);
		
		*deltae32=deltah32;
	}
	
	*f1 = t1-deltam12/n;
	*f2 = t3-deltam32/n;

	*q1 = sqrt(*f1* *f1+*f2* *f2);
	


 }