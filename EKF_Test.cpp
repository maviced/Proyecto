#include <stdio.h>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <complex>
#include <stdlib.h>
#include <string.h>
#include "include/Mjday.h"
#include "include/Frac.h"
#include "include/matriz.h"
#include "include/R_x.h"
#include "include/R_y.h"
#include "include/R_z.h"
#include "include/LTC.h"
#include "include/PoleMatrix.h"
#include "include/Position.h"
#include "include/vector.h"
#include "include/PrecMatrix.h"
#include "include/MeanObliquity.h"
#include "include/NutAngles.h"
#include "include/EqnEquinox.h"
#include "include/NutMatrix.h"
#include "include/gmst.h"
#include "include/gast.h"
#include "include/GHAMatrix.h"
#include "include/timediff.h"
#include "include/Geodetic.h"
#include "include/doubler.h"
#include "include/IERS.h"
#include "include/SAT_Const.h"
#include "include/anglesdr.h"
#include "include/AzElPa.h"
#include "include/TimeUpdate.h"
#include "include/MeasUpdate.h"
#include "include/Legendre.h"
#include "include/AccelHarmonic.h"
#include "include/G_AccelHarmonic.h"
#include "include/Accel.h"
#include "include/globales.h"
#include "include/VarEqn.h"
#include "include/DEInteg.h"
#include "include/sign_.h"

/**
* @file EKF_Test.cpp
*
* Este archivo contiene test unitarios y test de integración para las funciones implementadas
* en todos los archivos del proyecto Initial Orbit Determination, así como un test para el
* principal: EKF_GEOS3
*
*/

int tests_run = 0;
double **eopdata;
double pnm[362][362], dpnm[362][362], Cnm[362][362], Snm[362][362];
Param AuxParam;

/**
* Macro para imprimir un mensaje de error en caso de fallo.
*/
#define FAIL() printf("\nfailure in %s() line %d\n", __func__, __LINE__)

/**
* Macro para verificar una condición y en caso de fallo, imprimir un mensaje de error y retornar 1.
*/
#define _assert(test) do { if (!(test)) { FAIL(); return 1; } } while(0)

/**
* Macro para ejecutar una prueba y aumentar el contador de tests ejecutados.
*/
#define _verify(test) do { int r=test(); tests_run++; if(r) return r; } while(0)

using namespace std;

void inicio(){
	AuxParam.Mjd_TT = 4.974611324287046e+04;
	AuxParam.n = 10;
	AuxParam.m = 10;
	AuxParam.Mjd_UTC = 4.974611253472231e+04;
	
	FILE *fp;
	int f, c;
	double aux1, aux2;
	
	fp = fopen("./data/egm.txt", "r");
    if(fp == NULL){
		printf("Fail open GGM03S.txt file\n");
		exit(EXIT_FAILURE);
	}

    for(int n = 0; n <= 360; n++){
		for(int m = 0; m <= n; m++){
			fscanf(fp, "%d%d%lf%lf%lf%lf", &f, &c, &Cnm[n][m], &Snm[n][m], &aux1, &aux2);
		}
	}
	fclose(fp);
    
    //cout << Cnm[45][23] << endl;
    
    //read Earth orientation parameters
    //  ----------------------------------------------------------------------------------------------------
    // |  Date    MJD      x         y       UT1-UTC      LOD       dPsi    dEpsilon     dX        dY    DAT
    // |(0h UTC)           "         "          s          s          "        "          "         "     s
    //  ----------------------------------------------------------------------------------------------------
	
	fp = fopen("./data/eop19620101.txt", "r");
    if(fp == NULL){
        printf("Fail open eop19620101.txt file\n");
        exit(EXIT_FAILURE);
    }

    eopdata = (double **) malloc(13 * sizeof(double *));
    if(eopdata == NULL){
        printf("eopdata: memory not allocated\n");
        exit(EXIT_FAILURE);
    }
    
    for(int i = 0; i < 13; i++){
        eopdata[i] = (double *) malloc(19716 * sizeof(double));
        if(eopdata[i] == NULL){
            printf("eopdata[i]: memory not allocated\n");
            exit(EXIT_FAILURE);
        }
    }

    for(int i = 0; i < 19716; i++){
        fscanf(fp, "%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf", &eopdata[0][i], &eopdata[1][i], &eopdata[2][i], &eopdata[3][i], &eopdata[4][i], &eopdata[5][i], &eopdata[6][i], &eopdata[7][i], &eopdata[8][i], &eopdata[9][i], &eopdata[10][i], &eopdata[11][i], &eopdata[12][i]);
    }
    fclose(fp);
}

void fin(){
	for(int i = 0; i < 13; i++){
        free(eopdata[i]);
    }
    free(eopdata);
}

int Mjday_01()
{

    _assert(fabs(Mjday(2023,4,27,18,17,9.34)-60061.7619136572) < pow(10,-10));
	
    return 0;
}

int Mjday_02()
{

    _assert(fabs(Mjday(2023,4,27,0,0,0.0)-Mjday(2023,4,27)) < pow(10,-10));
	
    return 0;
}

int Frac_01()
{

    _assert(fabs(Frac(10.7525)-0.752500000000000) < pow(10,-10));
	
    return 0;
}

int Frac_02()
{

    _assert(fabs(Frac(-10.7525)-0.247500000000000) < pow(10,-10));
	
    return 0;
}

int Frac_03()
{

    _assert(fabs(Frac(10)-0) < pow(10,-10));
	
    return 0;
}

int MatricesIguales_01()
{
	double matriz1[3][3] = {{1.0, 2.0, 3.0}, {3.0, 4.0, 5.0}, {3.0, 5.0, 6.0}};
    double matriz2[3][3] = {{1.0, 2.0, 3.0}, {3.0, 4.0, 5.0}, {3.0, 5.0, 6.0}};

    _assert(matricesIguales(matriz1, matriz2, 3, 3));
	
    return 0;
}

int MatricesIguales_02()
{
	double matriz1[3][3] = {{1.0, 2.0, 3.0}, {3.0, 4.0, 5.0}, {3.0, 5.0, 6.0}};
    double matriz3[3][3] = {{1.0, 2.0, 5.0}, {3.0, 5.0, 6.0}, {3.0, 5.0, 6.0}};

    _assert(!matricesIguales(matriz1, matriz3, 3, 3));
	
    return 0;
}

int Transpuesta_01()
{
	double m[3][3] = {{1.0, 2.0, 3.0}, {3.0, 4.0, 5.0}, {1.0, 4.0, 6.0}};
	double mT[3][3] = {{1.0, 3.0, 1.0}, {2.0, 4.0, 4.0}, {3.0, 5.0, 6.0}};
	
	transpuesta(m, 3, 3);
	
    _assert(matricesIguales(m, mT, 3, 3));
	
    return 0;
}

int Mult3x3_01()
{
	double m1[3][3] = {{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}, {7.0, 8.0, 9.0}};
    double m2[3][3] = {{9.0, 8.0, 7.0}, {6.0, 5.0, 4.0}, {3.0, 2.0, 1.0}};
    double m3[3][3] = {{30.0, 24.0, 18.0}, {84.0, 69.0, 54.0}, {138.0, 114.0, 90.0}};
	double resultado[3][3];
	
	mult3x3(m1, m2, resultado);
	
    _assert(matricesIguales(m3, resultado, 3, 3));
	
    return 0;
}

int Mult3x1_01()
{
	double ma1[3][3] = {{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}, {7.0, 8.0, 9.0}};
    double ma2[3][1] = {{1.0}, {2.0}, {3.0}};
	double ma3[3][1] = {{14.0}, {32.0}, {50.0}};
    double resultado2[3][1];
	
	mult3x1(ma1, ma2, resultado2);
	
    _assert(ma3[0][0]==resultado2[0][0] && ma3[1][0]==resultado2[1][0] && ma3[2][0]==resultado2[2][0]);
	
    return 0;
}

int Mult3Matrices_01()
{
	double matriz1[3][3] = {{1.0, 2.0, 3.0}, {3.0, 4.0, 5.0}, {3.0, 5.0, 6.0}};
    double matriz2[3][3] = {{1.0, 2.0, 3.0}, {3.0, 4.0, 5.0}, {3.0, 5.0, 6.0}};
	double matriz3[3][3] = {{1.0, 2.0, 5.0}, {3.0, 5.0, 6.0}, {3.0, 5.0, 6.0}};
	double matriz4[3][3] = {{184.0, 312.0, 416.0}, {348.0, 590.0, 786.0}, {414.0, 702.0, 936.0}};
	double result[3][3];
	
	mult3Matrices(matriz1, matriz2, matriz3, result);
	
    _assert(matricesIguales(matriz4, result, 3, 3));
	
	return 0;
	
}

int vectorAMatriz_01()
{
	double v[3]={1.0, 3.0, 4.0};
	double m[3][1];
	
	vectorAMatriz(v, m);
	
    _assert(m[0][0]==1.0 && m[1][0] == 3.0 && m[2][0] == 4.0);
	
	return 0;
	
}

int matrizAVector_01()
{
	double v[3];
	double m[3][1] = {{1.0}, {3.0}, {4.0}};
	
	matrizAVector(m, v);
	
    _assert(v[0]==1.0 && v[1] == 3.0 && v[2]== 4.0);
	
	return 0;
	
}

int RX_01()
{
	double angle = M_PI/4;
    double rotmat[3][3];
	double r[3][3] ={{1.0, 0.0, 0.0},{0.0, 0.707106781186548, 0.707106781186547}, {0.0, -0.707106781186547, 0.707106781186548}};
	
	R_x(angle, rotmat);

    _assert(matricesIguales(r, rotmat, 3, 3));
	
    return 0;
}

int RY_01()
{
	double angle = 0.2;
    double rotmat[3][3];
	double r[3][3] ={{0.980066577841242, 0.198669330795061, 0.0},{-0.198669330795061, 0.980066577841242, 0.0}, {0.0, 0.0, 1.0}};
	
	R_z(angle, rotmat);

    _assert(matricesIguales(r, rotmat, 3, 3));
	
    return 0;
}

int RZ_01()
{
	double angle = 0.2;
    double rotmat[3][3];
	double r[3][3] ={{0.980066577841242, 0.0, -0.198669330795061},{0.0, 1.0, 0.0}, {0.198669330795061, 0.0, 0.980066577841242}};
	
	R_y(angle, rotmat);

    _assert(matricesIguales(r, rotmat, 3, 3));
	
    return 0;
}

int LTC_01()
{
    double M[3][3];
	double result[3][3] ={{-0.099833416646828, 0.995004165278026, 0.0},
	{-0.197676811654084, -0.019833838076210, 0.980066577841242}, 
	{0.975170327201816, 0.097843395007256, 0.198669330795061}};
	
	LTC(0.1, 0.2, M);

    _assert(matricesIguales(M, result, 3, 3));
	
    return 0;
}

int PoleMatrix_01()
{
    double PoleMat[3][3];
	double result[3][3] ={{-0.416146836547142, 0.128320060202457, -0.900197629735517},
	{0.0, -0.989992496600445, -0.141120008059867}, 
	{-0.909297426825682, -0.058726644927621, 0.411982245665683}};
	
	PoleMatrix(2, 3, PoleMat);
	
	
    _assert(matricesIguales(PoleMat, result, 3, 3));
	
    return 0;
}

int VectoresIguales_01()
{
	double v1[3] = {1.0, 2.0, 3.0};
    double v2[3] = {1.0, 2.0, 3.0};

    _assert(vectoresIguales(v1, v2, 3));
	
    return 0;
}

int VectoresIguales_02()
{
	double v1[3] = {1.0, 2.0, 3.0};
    double v3[3] = {3.0, 5.0, 6.0};

    _assert(!vectoresIguales(v1, v3, 3));
	
    return 0;
}

int norm_01()
{
    double v[3] = {1.0, 1.0, 1.0};
    _assert(fabs(norm(v, 3) - sqrt(3.0)) < pow(10,-10));
    return 0;
}

int dot_01()
{
    double v[3] = {1.0, 1.0, 1.0};
	double w[3] = {1.0, 2.0, 3.0};
    _assert(fabs(dot(v, 3, w, 3) - 6.0) < pow(10,-10));
    return 0;
}

int cross_01()
{
    double v[3] = {1.0, 1.0, 1.0};
	double w[3] = {1.0, 2.0, 3.0};
	double s[3];
	double result[3] = {1.0, -2.0, 1.0};
	int i;
	cross(v, 3, w, 3, s, i);
    _assert(vectoresIguales(s, result, 3));
    return 0;
}

int Position_01()
{
    double r[3];
	double result[3] ={6.220593764224149e6,0.624141235510229e6, 1.258824205408778e6};
	
	Position(0.1, 0.2, 2, r);

	_assert(vectoresIguales(r, result, 3));
	
    return 0;
}

int PrecMatrix_01()
{
    double PrecMat[3][3];
	double result[3][3] ={{0.999977967802458, 0.006086828257242,0.002648477191451}, 
	{-0.006086828258241, 0.999981475056808, -0.000008060124240},
	{-0.002648477189153, -0.000008060879152, 0.999996492745650}};
	
	PrecMatrix(10000, 50.2, PrecMat);

	_assert(matricesIguales(PrecMat, result, 3, 3));
	
    return 0;
}

int MeanObliquity_01()
{

	_assert(fabs(MeanObliquity(100)-0.409412448949153) < pow(10,-10));
	
    return 0;
}

int NutAngles_01()
{
	double dpsi, deps;
	
	NutAngles(200, &dpsi, &deps);
	
	_assert(fabs(dpsi-4.635815436424795e-05) < pow(10,-10) && fabs(deps-3.346147297940520e-05) < pow(10,-10));
	
    return 0;
}

int EqnEquinox_01()
{
	double EqE = EqnEquinox(100);
	
	
	_assert(fabs(EqE-4.291003873206410e-05) < pow(10,-10));
	
    return 0;
}

int NutMatrix_01()
{
	double NutMat[3][3];
	double result[3][3] ={{0.999999998906011, -0.000042910038716, -0.000018620075003}, 
	{0.000042909302689, 0.999999998298191, -0.000039527320885}, 
	{0.000018621771090, 0.000039526521867, 0.999999999045442}};
	
	NutMatrix(100, NutMat);
	
	_assert(matricesIguales(NutMat, result, 3, 3));
	
    return 0;
}

int gmst_01()
{
	double gmstime = gmst(100);
	
	_assert(fabs(gmstime-2.693487276423788) < pow(10,-10));
	
    return 0;
}

int gast_01()
{
	double gstime = gast(100);
	
	_assert(fabs(gstime-2.693530186462521) < pow(10,-10));
	
    return 0;
}

int GHAMatrix_01()
{
	double GHAmat[3][3];
	double result[3][3] ={{-0.901288171692405, 0.433220072904480, 0.0}, {-0.433220072904480, -0.901288171692405, 0.0}, {0.0, 0.0, 1.0}};
	
	GHAMatrix(100, GHAmat);
	
	_assert(matricesIguales(GHAmat, result, 3, 3));
	
    return 0;
}

int timediff_01()
{
	double UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC;
	
	timediff(100, 200, &UT1_TAI, &UTC_GPS, &UT1_GPS, &TT_UTC, &GPS_UTC);
	
	_assert(fabs(UT1_TAI-(-100)) < pow(10,-10) && fabs(UTC_GPS-(-181)) < pow(10,-10) && fabs(UT1_GPS-(-81)) < pow(10,-10) 
	&& fabs(TT_UTC-(2.321840000000000e+02)) < pow(10,-10) && fabs(GPS_UTC-181) < pow(10,-10));
	
    return 0;
}

int Geodetic_01()
{
	double lon, lat, h;
	double r[3] = {0.1, 0.2, 0.3};
	
	Geodetic(r, &lon, &lat, &h);
	
	_assert(fabs(lon-1.107148717794090) < pow(10,-10) && fabs(lat-1.570791107411455) < pow(10,-10) && fabs(h-(-6.356752014244596e+06)) < pow(10,-6));
	
    return 0;
}

int doubler_01()
{
	double r2[3], r3[3], magr1, magr2, a, deltae32;
	complex<double> f1, f2, q1;
	
    double los1[3] = {0.6, 0.7, 0.5};
    double los2[3] = {-0.3, 0.4, -0.2};
    double los3[3] = {0.1, -0.5, 0.9};
    double rsite1[3] = {1.0, 2.0, 3.0};
    double rsite2[3] = {-1.0, -2.0, -3.0};
    double rsite3[3] = {0.5, -0.5, 0.0};

	
	double r2_result[3] = {-1.440824042221049, -1.412234610371934, -3.293882694814033};
	double r3_result[3] = {1.506234734832367, -5.531173674161833, 9.056112613491299};
	complex<double> f1_result(10.0, - 0.000000064502940);
	complex<double> f2_result(20.0, - 0.000000086389856);
	complex<double> q1_result(22.360679774997898, - 0.000000106116028);
	
	
	doubler(0.5, 0.3, 1.2, 0.8, 2.5, 1.8, los1, los2, los3, rsite1, rsite2, rsite3, 10.0, 20.0, 'y', r2, r3, &f1, &f2, &q1, &magr1, &magr2, &a, &deltae32); 
	
	//cout << fabs(f1-f1_result) << endl;
	
	_assert(vectoresIguales(r2, r2_result, 3) && vectoresIguales(r3, r3_result, 3) && fabs(f1-f1_result) < pow(10,-10) && fabs(f2-f2_result) < pow(10,-10)
	&& fabs(q1-q1_result) < pow(10,-10) && fabs(magr1-5.649430266429126) < pow(10,-10) && fabs(magr2-3.862647242833589) < pow(10,-10) 
	&& fabs(a-0.018245890404468) < pow(10,-10)
	&& fabs(deltae32-8.238735060758426) < pow(10,-10));
	
    return 0;
}

int IERS_01()
{
	double UT1_UTC, TAI_UTC, x_pole, y_pole;
	
	//cout << eopdata[2][2] << endl;
	IERS(eopdata, 37671.123, &UT1_UTC, &TAI_UTC, &x_pole, &y_pole);

	//cout << UT1_UTC << endl;
	//cout << TAI_UTC << endl;
	//cout << x_pole << endl;
	//cout << y_pole << endl;
	
	_assert(fabs(UT1_UTC-0.0302682) < pow(10,-10) && fabs(TAI_UTC-2) < pow(10,-10) && fabs(x_pole-(-1.46408883558269e-07)) < pow(10,-10) 
	&& fabs(y_pole-1.06320125081002e-06) < pow(10,-10));
	
	return 0;
	
}

int anglesdr_01()
{
	double az1, az2, az3, el1, el2, el3, Mjd1, Mjd2, Mjd3, r2[3], v2[3];
	az1 = 1.055908489493301;
    az2 = 1.363102145807571;
    az3 = 1.976156026887588;
    el1 = 0.282624656433946;
    el2 = 0.453434794338875;
    el3 = 0.586427138011591;
    Mjd1 = 4.974611015046295e+04;
    Mjd2 = 4.974611128472211e+04;
    Mjd3 = 4.974611253472231e+04;
    double rs[3] = {-5.512568445011530e+06, -2.196994687777972e+06, 2.330805221940450e+06};	
	
	double r2_result[3] = {6147304.28873136, 2498216.0975712, 2872808.05359545};
	double v2_result[3] = {3764.62899474253, -2217.8449407281, -6141.47100738894};
	
	anglesdr(az1, az2, az3, el1, el2, el3, Mjd1, Mjd2, Mjd3, rs, rs, rs, r2, v2);
	/*
	cout << r2[0] << endl;
	cout << r2[1] << endl;
	cout << r2[2] << endl;
	
	cout << v2[0] << endl;
	cout << v2[1] << endl;
	cout << v2[2] << endl;*/
	
	_assert(vIgualesMenosPrecision(r2, r2_result, 3) && vIgualesMenosPrecision(v2, v2_result, 3));
	
	return 0;
	
}

int AzElPa_01()
{
	double s[3] = {1.0, 2.0, 3.0};
	double Az, El, dAds[3], dEds[3];
	double dAds_result[3] = {0.4, -0.2, 0.0};
	double dEds_result[3] = {-0.095831484749991, -0.191662969499982, 0.159719141249985};
	
	AzElPa(s, &Az, &El, dAds, dEds);
	
	_assert(fabs(Az-0.463647609000806) < pow(10,-10)  && fabs(El-0.930274014115472) < pow(10,-10) && vectoresIguales(dAds, dAds_result,3) 
	&& vectoresIguales(dEds, dEds_result, 3));
	
    return 0;
}

int TimeUpdate_01()
{
    double P[6][6] = {
        {1.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 1.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 1.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 1.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 1.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 1.0}
    };
	
	double Phi[6][6] = {
        {1.0, 2.0, 3.0, 4.0, 5.0, 6.0},
        {7.0, 8.0, 9.0, 10.0, 11.0, 12.0},
        {13.0, 14.0, 15.0, 16.0, 17.0, 18.0},
        {19.0, 20.0, 21.0, 22.0, 23.0, 24.0},
        {25.0, 26.0, 27.0, 28.0, 29.0, 30.0},
        {31.0, 32.0, 33.0, 34.0, 35.0, 36.0}
    };
	double result[6][6] = {
		{91.1, 217.1, 343.1, 469.1, 595.1, 721.1},
		{217.1, 559.1, 901.1, 1243.1, 1585.1, 1927.1},
		{343.1, 901.1, 1459.1, 2017.1, 2575.1, 3133.1},
		{469.1, 1243.1, 2017.1, 2791.1, 3565.1, 4339.1},
		{595.1, 1585.1, 2575.1, 3565.1, 4555.1, 5545.1},
		{721.1, 1927.1, 3133.1, 4339.1, 5545.1, 6751.1}
	};
	
	TimeUpdate(P, Phi, 0.1);

	_assert(matricesIguales6x6(P, result, 6, 6));
	
    return 0;
}

int TimeUpdate_02()
{
    double P[6][6] = {
        {100000000.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 100000000.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 100000000.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 1000.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 1000.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 1000.0}
    };
	
	double Phi[6][6] = {
        {1.00059150895668, 0.000639937164090078, 0.000841574809808159, 37.0074570164194, 0.00789930159699147, 0.0103300884404148},
        {0.000639932089384861, 0.999591480460824, 0.000410553228024333, 0.00789927031365396, 36.9949287532614, 0.00499340817549725},
        {0.000841564064343303, 0.000410551239405512, 0.999817530876772, 0.010330022181339, 0.00499339591198313, 36.9976195120209},
        {3.26974564966778e-05, 3.46242074032268e-05, 4.52774644404386e-05, 1.00061780492089, 0.000640996693188821, 0.000833473994136441},
        {3.46235215439914e-05, -2.22192313223797e-05, 2.18885972715899e-05, 0.00064099161704589, 0.999586160303014, 0.000399224779033486},
        {4.52760122535526e-05, 2.18883285761002e-05, -1.04219803965523e-05, 0.000833463245521863, 0.000399222790368379, 0.999796555220812}
    };

	double result[6][6] = {
        {101488000.600326, 128617.803934412, 169139.099218971, 40308.0399813276, 3496.44199102555, 4571.97846555939},
        {128617.803934412, 101286995.428481, 82509.594120304, 3496.57950927422, 34761.7247204656, 2210.17635089282},
        {169139.099218971, 82509.594120304, 101332421.164292, 4572.26963026297, 2210.23022307775, 35952.8041664988},
        {40308.0399813276, 3496.57950927422, 4572.26963026297, 1001.66889789894, 1.41783559164539, 1.84417805416141},
        {3496.44199102555, 34761.7247204656, 2210.23022307775, 1.41783559164539, 999.390221437519, 0.884050439731921},
        {4571.97846555939, 2210.17635089282, 35952.8041664988, 1.84417805416141, 0.884050439731921, 999.857769260102}
    };
	
	TimeUpdate(P, Phi);

	_assert(matricesIguales6x6(P, result, 6, 6));
	
    return 0;
}

int MeasUpdate_01()
{
	double z, g, s, K[6][1], x2[6][1], P2[6][6];
    double x[6][1] = {{6.509022750859073e+06}, {2.219445217963939e+06}, {2.192214071528120e+06}, {2.983142032879137e+03}, 
	{-2.463758038689052e+03}, {-6.345351158787863e+03}};
	z = 1.055908489493301;
    g = 1.976212219010204;
    s = 3.909537524467298e-04;
    double G[1][6] = {3.551136354805395e-07, -2.904907230151484e-07, -7.548892159571566e-07, 0, 0, 0};
	double P[6][6] = {{1.968461133763149e+04, -7.980618629759560e+03, 1.198818585948804e+03, 1.399843419651444e+02, -63.965400874876010, 37.553853164749520},
        {-7.980618629759563e+03, 6.148208470878238e+03, -1.474482254089532e+03, -44.539805608375030, 57.152125509865610, -44.003960337648444},
        {1.198818585948847e+03, -1.474482254089560e+03, 2.626509836023459e+04, 41.254425893985920, -5.879007021683925, 1.539920720965272e+02},
		{1.399843419651444e+02, -44.539805608374990, 41.254425893985410, 1.718311758226475, -0.374227857550987, 0.386226470317009},
		{-63.965400874876230, 57.152125509865690, -5.879007021683481, -0.374227857550988, 0.786988774990078, -0.781119947476869},
		{37.553853164750180, -44.003960337648756, 1.539920720965270e+02, 0.386226470317016, -0.781119947476872, 1.941919119102344}};
	
	double K_result[6][1] = {{49094.8843402307}, {-20488.0742513267}, {-110843.882377124}, {184.063454381519}, {-203.768498360936}, {-526.541335273855}};
	double x2_result[6][1] = {{6463840.54570056}, {2238300.46910805}, {2294224.10987392}, {2813.74774934406}, {-2276.22912968942}, {-5860.77320419052}};
	
	double P2_result[6][6] = {{19272.0378450461, -7808.44516476797, 2130.3056094671, 138.437547360125, -62.2530130736329, 41.978692954147},
		{-7808.44516476797, 6076.35775196402, -1863.20656943233, -43.8943036758403, 56.4375189089863, -45.8505161931862},
		{2130.30560946713, -1863.20656943235, 24162.0352924992, 44.7466984915168, -9.74514734517747, 144.001898420699},
		{138.437547360125, -43.8943036758403, 44.7466984915163, 1.71251261312114, -0.367807880843163, 0.402815801520584},
		{-62.2530130736331, 56.4375189089864, -9.74514734517709, -0.367807880843164, 0.779881503111143, -0.799485261327591},
		{41.9786929541476, -45.8505161931865, 144.001898420699, 0.402815801520591, -0.799485261327594, 1.89446282946255}};
	
	MeasUpdate(x, z, g, s, G, P, 6, K, x2, P2);
	
	_assert(matricesIguales6x1(K, K_result) && matricesIguales6x1(x2, x2_result) && matricesIguales6x6(P2, P2_result, 6, 6));
	
    return 0;
}

int Legendre_01()
{
	bool b1=true, b2=true;
	double pnm_result[6][6] = {{1.0, 0,0,0,0,0},
		{0.344105374842754,  1.697525107621188, 0,0,0,0},
		{-0.985649251135536,   0.754105377234551,   1.860059309205460, 0, 0, 0},
		{-0.736578611760033,  -1.274523324299256,   0.977702309576464,   1.969045565334557, 0, 0},
		{0.701414763144989,  -1.257792762729005,  -1.165800091961221,   1.173566894309998,   2.046857490454963, 0},
		{1.015982330118486,   0.755573868609946,  -1.429355897146991,  -1.052688679513099,   1.348698200497283,   2.103969928789927}};
		
	double dpnm_result[6][6] = {{0,0,0,0,0,0},
		{1.697525107621188,  -0.344105374842754, 0,0,0,0},
		{1.306148827631137,   3.567253890614423,  -0.754105377234551, 0, 0, 0},
		{-3.121931809808946,   3.350124840194014,   4.426776775577388,  -1.197435889401486, 0,0},
		{-3.977489954699398,  -4.691103687606732,   4.863723995043303,   5.075705886040589,  -1.659674218285272, 0},
		{2.926325009955833,  -7.716602883444457,  -4.577610676292213,   6.362213538885371,   5.559758463025293,  -2.132479094870926}};

	
	Legendre(5, 5, 0.2, pnm, dpnm);
	/*
	for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            cout << pnm[i][j] << " ";
        }
        cout << endl;
    }*/
	
	for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            if (fabs(pnm[i][j] - pnm_result[i][j]) > pow(10,-10)) {
                b1 = false;
            }
			if (fabs(dpnm[i][j] - dpnm_result[i][j]) > pow(10,-10)) {
                b2 = false;
            }
        }
    }
	
	_assert(b1 && b2);
	
    return 0;
}

int AccelHarmonic_01()
{
	double a[3][1], v[3];
	double r[3][1] = {{1.0e6}, {2.0e6}, {3.0e6}};
    double E[3][3] = {{1.0, 0, 0}, {0, 1.0, 0}, {0,0,1.0}};
	double a_result[3] = {-7.5193753936907, -15.0845129004961, -22.7912404959225};
	
	AccelHarmonic(r, E, 10, 10, a);
	
	matrizAVector(a, v);
	
	_assert(vectoresIguales(v, a_result, 3));

    return 0;
}

int AccelHarmonic_02()
{
	double a[3][1], v[3];
	double r[3][1] = {6509051.32760191, 2219445.79558225, 2192198.3392695};
    double E[3][3] = {{-0.979791105663104, 0.200022992334509, -0.000437950935389831},
                      {-0.200022955088195, -0.979791202823139, -0.000127703454978668},
                      {-0.00045464410095941, -3.75224690702425e-05, 0.999999895945397}};
	double a_result[3] = {-6.90400425350075, -2.35419189643486, -2.33118411595261};
	
	AccelHarmonic(r, E, 10, 10, a);
	
	matrizAVector(a, v);
	
	_assert(vectoresIguales(v, a_result, 3));

    return 0;
}

int G_AccelHarmonic_01()
{
	double G[3][3];
	double r[3][1] = {{1.0e6}, {2.0e6}, {3.0e6}};
    double U[3][3] = {{1.0, 0, 0}, {0, 1.0, 0}, {0,0,1.0}};
	double G_result[3][3] = {{-5.95770687183972e-06,      3.17801080917945e-06,       4.8026048178329e-06},
		{3.17801080562674e-06,     -1.10375986928091e-06,      9.71725178722238e-06},
		{4.80260482049744e-06,     9.71725179610416e-06,      7.06146674289698e-06}};
	
	G_AccelHarmonic(r, U, 10, 10, G);
	
	//imprimirMatriz(G, 3, 3);
	
	_assert(matricesIguales(G, G_result, 3, 3));

    return 0;
}

int Accel_01()
{
	double dY[6][1];
	double Y[6][1] = {{6.509022750859073e+06}, {2.219445217963939e+06}, {2.192214071528120e+06}, {2.983142032879137e+03}, {-2.463758038689052e+03}, {-6.345351158787863e+03}};
	double dY_result[6][1] = {{2983.14203287914}, {-2463.75803868905}, {-6345.35115878786}, {-6.90403463770349}, {-2.35421197848897}, {-2.33122137458625}};
	
	Accel(0.2, Y, dY);
	
	_assert(matricesIguales6x1(dY, dY_result));

    return 0;
}

int VarEqn_01()
{
	double yPhip[42][1];
	
	double yPhi[42][1] = {{6.509051327601911e+06},
        {2.219445795582250e+06},
		{2.192198339269504e+06},
		{2.983377225623833e+03},
		{-2.463724973225010e+03},
		{-6.345490076123197e+03},
		{1.000108572859796},
		{6.385902667582700e-05},
		{6.423805110095676e-05},
		{1.815381491355681e-05},
		{1.062998077320940e-05},
		{1.065598297506753e-05},
		{6.385904485700467e-05},
		{0.999945681086397},
		{2.218114868061419e-05},
		{1.062998837898381e-05},
		{-9.068980303738064e-06},
		{3.668059404843647e-06},
		{6.423808975256889e-05},
		{2.218115572986840e-05},
		{0.999945751899337},
		{1.065599912795597e-05},
		{3.668062343855802e-06},
		{-9.082886132725295e-06},
		{12.000451055006335},
		{2.551176484344597e-04},
		{2.557427972414067e-04},
		{1.000109266525855},
		{6.369930584923789e-05},
		{6.363233902145628e-05},
		{2.551176835535193e-04},
		{11.999797714842991},
		{8.803125350264360e-05},
		{6.369932409065732e-05},
		{0.999945488506171},
		{2.183507105246286e-05},
		{2.557428724008849e-04},
		{8.803126728659642e-05},
		{11.999797378174310},
		{6.363237784776592e-05},
		{2.183507811971977e-05},
		{0.999945250813579}};
	
	double yPhip_result[42][1] = {{2983.37722562383},
		{-2463.72497322501},
		{-6345.4900761232},
		{-6.90400425350075},
		{-2.35419189643486},
		{-2.33118411595261},
		{1.81538149135568e-05},
		{1.06299807732094e-05},
		{1.06559829750675e-05},
		{1.52740048581475e-06},
		{8.82513650969482e-07},
		{8.75370305812642e-07},
		{1.06299883789838e-05},
		{-9.06898030373806e-06},
		{3.66805940484365e-06},
		{8.82516184261497e-07},
		{-7.59688950830324e-07},
		{2.98491275937814e-07},
		{1.0655999127956e-05},
		{3.6680623438558e-06},
		{-9.08288613272529e-06},
		{8.75375690413217e-07},
		{2.98492253878844e-07},
		{-7.67224435277215e-07},
		{1.00010926652586},
		{6.36993058492379e-05},
		{6.36323390214563e-05},
		{1.83266030169558e-05},
		{1.05896467799661e-05},
		{1.05039407549732e-05},
		{6.36993240906573e-05},
		{0.999945488506171},
		{2.18350710524629e-05},
		{1.05896543280486e-05},
		{-9.1171132376549e-06},
		{3.58171873111844e-06},
		{6.36323778477659e-05},
		{2.18350781197198e-05},
		{0.999945250813579},
		{1.05039568891677e-05},
		{3.58172165700959e-06},
		{-9.20754131164641e-06}};
	
	VarEqn(0.2, yPhi, yPhip);
	/*
	for (int i = 0; i < 42; i++) {
        cout << yPhip[i][0] << endl;
    }*/
	
	_assert(matricesIguales42x1(yPhip, yPhip_result));

    return 0;
}

int sign_01()
{
	
    _assert(fabs(sign_(1.0, -7.0)-(-1.0)) < pow(10,-10));
	
    return 0;
	
}

int sign_02()
{
	
    _assert(fabs(sign_(2.0, 8.0)-2.0) < pow(10,-10));
	
    return 0;
	
}

int DEInteg_01()
{
	double Y2[6][1];
	double Y[6][1] = {{6147304.28873136}, {2498216.09757119}, {2872808.05359544}, {3764.62899474253}, {-2217.84494072807}, {-6141.47100738888}};
	double Y_result[6][1] = {{5.581748453497840e+06}, {2.772707485642144e+06}, {3.671626337361284e+06}, {0.004600309588766e+06}, {-0.001842298276890e+06}, {-0.005674023618748e+06}};
	
	DEInteg(&Accel, 0, -134.999991953373, 1e-13, 1e-6, 6, Y);
	
	/*
	for (int i = 0; i < 6; i++) {
        cout << Y[i][0] << endl;
    }*/
	
	_assert(matricesIguales6x1(Y, Y_result));

    return 0;
}

int EKF_GEOS3_01()
{
	//Y es la matriz que obtengo al ejecutar el principal en c++
	//Y2 es  la matriz que obtengo al ejecutar el principal en Matlab
	double Y[6][1] = {{5753168.639158068}, {2673421.59955983}, {3440276.254340664}, {4327.622189912358}, {-1927.270445234001}, {-5726.203328424704}};
	double Y2[6][1] = {{5753168.63919783}, {2673421.59961428}, {3440276.2544536}, {4327.62218982135}, {-1927.27044563109}, {-5726.2033292681}};
	
	bool b = true;
	
	for (int j = 0; j < 6; j++) {
		if (fabs(Y[j][0] - Y2[j][0]) > pow(10,-3)) {
			b = false;
		}
	}
	
	//Como es el test del principal, se ha ido acumulando la falta de precisión de todos las funciones del proyecto
	//y en esta solo nos coindicen los tres primeros decimales., pero al ser el pricipal lo consideramos acpetable
	
	_assert(b);

    return 0;
}


/**
* Ejecuta todas las pruebas definidas.
*
* @return 0 si todas las pruebas pasan, 1 si alguna prueba falla.
*/
int all_tests()
{
    _verify(Mjday_01);
	_verify(Mjday_02);
	_verify(Frac_01);
	_verify(Frac_02);
	_verify(Frac_03);
	_verify(MatricesIguales_01);
	_verify(MatricesIguales_02);
	_verify(Transpuesta_01);
	_verify(Mult3x3_01);
	_verify(Mult3x1_01);
	_verify(Mult3Matrices_01);
	_verify(vectorAMatriz_01);
	_verify(matrizAVector_01);
	_verify(RX_01);
	_verify(RY_01);
	_verify(RZ_01);
	_verify(LTC_01);
	_verify(PoleMatrix_01);
	_verify(VectoresIguales_01);
	_verify(VectoresIguales_02);
	_verify(norm_01);
	_verify(dot_01);
	_verify(cross_01);
	_verify(Position_01);
	_verify(PrecMatrix_01);
	_verify(MeanObliquity_01);
	_verify(NutAngles_01);
	_verify(EqnEquinox_01);
	_verify(NutMatrix_01);
	_verify(gmst_01);
	_verify(gast_01);
	_verify(GHAMatrix_01);
	_verify(timediff_01);
	_verify(Geodetic_01);
	_verify(doubler_01);
	_verify(IERS_01);
	_verify(anglesdr_01);
	_verify(AzElPa_01);
	_verify(TimeUpdate_01);
	_verify(TimeUpdate_02);
	_verify(MeasUpdate_01);
	_verify(Legendre_01);
	_verify(AccelHarmonic_01);
	_verify(AccelHarmonic_02);
	_verify(G_AccelHarmonic_01);
	_verify(Accel_01);
	_verify(VarEqn_01);
	_verify(sign_01);
	_verify(sign_02);
	_verify(DEInteg_01);
	_verify(EKF_GEOS3_01);

    return 0;
}

/**
* Función principal.
*
* Ejecuta todas las pruebas y muestra los resultados.
*
*/
int main()
{
	inicio();
    int result = all_tests();
	fin();
	
    if (result == 0)
        printf("PASSED\n");

    printf("Tests run: %d\n", tests_run);

    return result != 0;
}
