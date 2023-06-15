//$Header$
//------------------------------------------------------------------------------
// matriz
//------------------------------------------------------------------------------
// EKF_GEOS3: Initial Orbit Determination using Double-R-Iteration and Extended 
// Kalman Filter methods
//
//
// Last modified:   2015/08/12   M. Mahooti
//
/**
*
* @file matriz.h
*
* Provides basic operations to manipulate static matrixs of different sizes.
*
*/
//------------------------------------------------------------------------------

#ifndef _IMPRIMIRMATRIZ_
#define _IMPRIMIRMATRIZ_

void imprimirMatriz(double matriz[][3], int filas, int columnas);

#endif


#ifndef _MATRICESIGUALES_
#define _MATRICESIGUALES_


bool matricesIguales(double matriz1[][3], double matriz2[][3], int filas, int columnas);

#endif


#ifndef _TRANSPUESTA_
#define _TRANSPUESTA_

void transpuesta(double matriz[][3], int filas, int columnas);

#endif


#ifndef _MULT3X3_
#define _MULT3X3_

void mult3x3(double matriz1[3][3], double matriz2[3][3], double resultado[3][3]);

#endif


#ifndef _MULT3X1_
#define _MULT3X1_

void mult3x1(double matriz1[3][3], double matriz2[3][1], double resultado[3][1]);

#endif

#ifndef _MULT3MATRICES_
#define _MULT3MATRICES_

void mult3Matrices(double matriz1[3][3], double matriz2[3][3], double matriz3[3][3], double resultado[3][3]);

#endif

#ifndef _vectorAMatriz_
#define _vectorAMatriz_

void vectorAMatriz(double v[3], double m[3][1]);

#endif

#ifndef _matrizAVector_
#define _matrizAVector_

void matrizAVector(double m[3][1], double v[3]);

#endif

#ifndef _MULT6X6_
#define _MULT6X6_

void mult6x6(double matriz1[6][6], double matriz2[6][6], double resultado[6][6]);

#endif

#ifndef _MULT3MATRICES6x6_
#define _MULT3MATRICES6x6_

void mult3Matrices6x6(double matriz1[6][6], double matriz2[6][6], double matriz3[6][6], double resultado[6][6]);

#endif

#ifndef _TRANSPUESTA6x6_
#define _TRANSPUESTA6x6_

void transpuesta6x6(double matriz[6][6], int filas, int columnas);

#endif

#ifndef _copiarMatriz6x6_
#define _copiarMatriz6x6_

void copiarMatriz6x6(double matriz1[6][6], double matriz2[6][6]);

#endif

#ifndef _MATRICESIGUALES6x6_
#define _MATRICESIGUALES6x6_


bool matricesIguales6x6(double matriz1[6][6], double matriz2[6][6], int filas, int columnas);

#endif

#ifndef _multMatrices1x6_6x6_
#define _multMatrices1x6_6x6_


void multMatrices1x6_6x6(double matriz1[1][6], double matriz2[6][6], double result[1][6]);

#endif

#ifndef _transpuesta1x6_
#define _transpuesta1x6_

void transpuesta1x6(double matriz[1][6], double result[6][1]);

#endif

#ifndef _multMatrices1x6_6x1_
#define _multMatrices1x6_6x1_

double multMatrices1x6_6x1(double matriz1[1][6], double matriz2[6][1]);

#endif

#ifndef _multMatrices6x6_6x1_
#define _multMatrices6x6_6x1_

void multMatrices6x6_6x1(double matriz1[6][6], double matriz2[6][1], double result[6][1]);

#endif

#ifndef _eye6x6_
#define _eye6x6_

void eye6x6(double matriz[6][6]);

#endif

#ifndef _multMatrices6x1_1x6_
#define _multMatrices6x1_1x6_

void multMatrices6x1_1x6(double matriz1[6][1], double matriz2[1][6], double result[6][6]);

#endif

#ifndef _multMatrices6x6_
#define _multMatrices6x6_

void multMatrices6x6(double matriz1[6][6], double matriz2[6][6], double result[6][6]);

#endif

#ifndef _MATRICESIGUALES6x1_
#define _MATRICESIGUALES6x1_

bool matricesIguales6x1(double matriz1[6][1], double matriz2[6][1]);

#endif

#ifndef _MATRICESIGUALES42x1_
#define _MATRICESIGUALES42x1_

bool matricesIguales42x1(double matriz1[42][1], double matriz2[42][1]);

#endif

#ifndef _mult1x3_
#define _mult1x3_

void mult1x3(double matriz1[1][3], double matriz2[3][3], double resultado[1][3]);

#endif

#ifndef _vectorAMatriz2_
#define _vectorAMatriz2_

void vectorAMatriz2(double v[3], double m[1][3]);

#endif