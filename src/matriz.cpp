#include <iostream>
#include <cmath>
#include "../include/matriz.h"

/**
*
* @file matriz.cpp
*
* Provides basic operations to manipulate static matrixs of different sizes.
*
*/

using namespace std;

//------------------------------------------------------------------------------
// void imprimirMatriz(double matriz[][3], int filas, int columnas)
//------------------------------------------------------------------------------
/**
*
* This function prints the elements of a matrix in row-major order.
*
* @param matriz The matrix to be printed.
* @param filas The number of rows in the matrix.
* @param columnas The number of columns in the matrix.
*/ 
//------------------------------------------------------------------------------

void imprimirMatriz(double matriz[][3], int filas, int columnas) {
    for (int i = 0; i < filas; i++) {
        for (int j = 0; j < columnas; j++) {
            cout << matriz[i][j] << " ";
        }
        cout << endl;
    }
}

//------------------------------------------------------------------------------
// bool matricesIguales(double matriz1[][3], double matriz2[][3], int filas, int columnas)
//------------------------------------------------------------------------------
/**
*
* This function compares the elements of the two matrices and returns true if all the 
* elements are equal within a tolerance of 10^-10, and false otherwise.
*
* @param matriz1 The first matrix to compare.
* @param matriz2 The second matrix to compare.
* @param filas The number of rows in the matrices.
* @param columnas The number of columns in the matrices.
* @return True if the matrices are equal, false otherwise.
*/ 
//------------------------------------------------------------------------------

bool matricesIguales(double matriz1[][3], double matriz2[][3], int filas, int columnas) {
    
    for (int i = 0; i < filas; i++) {
        for (int j = 0; j < columnas; j++) {
            if (fabs(matriz1[i][j] - matriz2[i][j]) > pow(10,-10)) {
                return false;
            }
        }
    }
    return true;
}

//------------------------------------------------------------------------------
// void transpuesta(double matriz[][3], int filas, int columnas)
//------------------------------------------------------------------------------
/**
*
* This function transposes a matrix in-place by swapping elements across the main diagonal.
*
* @param matriz The matrix to be transposed.
* @param filas The number of rows in the matrix.
* @param columnas The number of columns in the matrix.
*/ 
//------------------------------------------------------------------------------

void transpuesta(double matriz[][3], int filas, int columnas) {
    double aux;
    for (int i = 0; i < filas; i++) {
        for (int j = i+1; j < columnas; j++) {
            aux = matriz[i][j];
            matriz[i][j] = matriz[j][i];
            matriz[j][i] = aux;
        }
    }
}


//------------------------------------------------------------------------------
// void mult3x3(double matriz1[3][3], double matriz2[3][3], double resultado[3][3])
//------------------------------------------------------------------------------
/**
*
* This function performs matrix multiplication between two 3x3 matrices and stores the result in a third 3x3 matrix.
*
* @param matriz1 The first matrix to be multiplied.
* @param matriz2 The second matrix to be multiplied.
* @param resultado The matrix to store the result of the multiplication.
*/ 
//------------------------------------------------------------------------------

void mult3x3(double matriz1[3][3], double matriz2[3][3], double resultado[3][3]) {
    for(int i = 0; i < 3; i++) {
        for(int j = 0; j < 3; j++) {
            resultado[i][j] = 0;
            for(int k = 0; k < 3; k++) {
                resultado[i][j] += matriz1[i][k] * matriz2[k][j];
            }
        }
    }
}


//------------------------------------------------------------------------------
// void mult3x1(double matriz1[3][3], double matriz2[3][1], double resultado[3][1])
//------------------------------------------------------------------------------
/**
*
* This function performs matrix multiplication between a 3x3 matrix and a 3x1 matrix and stores the result in a 3x1 matrix.
*
* @param matriz1 The 3x3 matrix to be multiplied.
* @param matriz2 The 3x1 matrix to be multiplied.
* @param resultado The matrix to store the result of the multiplication.
*/ 
//------------------------------------------------------------------------------

void mult3x1(double matriz1[3][3], double matriz2[3][1], double resultado[3][1]) {
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 1; j++) {
			resultado[i][j] = 0;
			for (int k = 0; k < 3; k++) {
				resultado[i][j] += matriz1[i][k] * matriz2[k][j];
			}
		}
        
    }
}


//------------------------------------------------------------------------------
// void mult3Matrices(double matriz1[3][3], double matriz2[3][3], double matriz3[3][3], double resultado[3][3])
//------------------------------------------------------------------------------
/**
*
* This function performs matrix multiplication between three 3x3 matrices and stores the result in a 3x3 matrix.
*
* @param matriz1 The first matrix to be multiplied.
* @param matriz2 The second matrix to be multiplied.
* @param matriz3 The third matrix to be multiplied.
* @param resultado The matrix to store the result of the multiplication.
*/ 
//------------------------------------------------------------------------------

void mult3Matrices(double matriz1[3][3], double matriz2[3][3], double matriz3[3][3], double resultado[3][3]) {
    double result[3][3];
	mult3x3(matriz1, matriz2, result);
	mult3x3(result, matriz3, resultado);
}

//------------------------------------------------------------------------------
// void vectorAMatriz(double v[3], double m[3][1])
//------------------------------------------------------------------------------
/**
*
* This function takes a 3-dimensional vector and converts it into a 3x1 matrix by assigning each vector element to the corresponding matrix element.
*
* @param v The 3-dimensional vector.
* @param m The 3x1 matrix to store the converted vector.
*/ 
//------------------------------------------------------------------------------

void vectorAMatriz(double v[3], double m[3][1]){
	for(int i=0; i<3; i++){
		m[i][0] = v[i];
	}
}

//------------------------------------------------------------------------------
// void matrizAVector(double m[3][1], double v[3])
//------------------------------------------------------------------------------
/**
*
* This function takes a 3x1 matrix and converts it into a 3-dimensional vector by assigning each matrix element to the corresponding vector element.
*
* @param m The 3x1 matrix to be converted.
* @param v The 3-dimensional vector to store the converted matrix.
*/ 
//------------------------------------------------------------------------------

void matrizAVector(double m[3][1], double v[3]){
	for (int i=0;i<3; i++){
		v[i] = m[i][0];
	}
}

//------------------------------------------------------------------------------
// void mult6x6(double matriz1[6][6], double matriz2[6][6], double resultado[6][6])
//------------------------------------------------------------------------------
/**
*
* This function performs matrix multiplication between two 6x6 matrices and stores the result in a 6x6 matrix.
*
* @param matriz1 The first matrix to be multiplied.
* @param matriz2 The second matrix to be multiplied.
* @param resultado The matrix to store the result of the multiplication.
*/ 
//------------------------------------------------------------------------------

void mult6x6(double matriz1[6][6], double matriz2[6][6], double resultado[6][6]) {
    for(int i = 0; i < 6; i++) {
        for(int j = 0; j < 6; j++) {
            resultado[i][j] = 0;
            for(int k = 0; k < 6; k++) {
                resultado[i][j] += matriz1[i][k] * matriz2[k][j];
            }
        }
    }
}

//------------------------------------------------------------------------------
// void mult3Matrices6x6(double matriz1[6][6], double matriz2[6][6], double matriz3[6][6], double resultado[6][6])
//------------------------------------------------------------------------------
/**
*
* This function performs matrix multiplication between three 6x6 matrices and stores the result in a 6x6 matrix.
*
* @param matriz1 The first matrix to be multiplied.
* @param matriz2 The second matrix to be multiplied.
* @param matriz3 The third matrix to be multiplied.
* @param resultado The matrix to store the result of the multiplication.
*/ 
//------------------------------------------------------------------------------

void mult3Matrices6x6(double matriz1[6][6], double matriz2[6][6], double matriz3[6][6], double resultado[6][6]) {
    double result[6][6];
	mult6x6(matriz1, matriz2, result);
	mult6x6(result, matriz3, resultado);
}

//------------------------------------------------------------------------------
// void transpuesta6x6(double matriz[6][6], int filas, int columnas)
//------------------------------------------------------------------------------
/**
*
* This function transposes a 6x6 matrix by swapping elements along the main diagonal.
*
* @param matriz The matrix to be transposed.
* @param filas The number of rows in the matrix.
* @param columnas The number of columns in the matrix.
*/ 
//------------------------------------------------------------------------------

void transpuesta6x6(double matriz[6][6], int filas, int columnas) {
    double aux;
    for (int i = 0; i < filas; i++) {
        for (int j = i+1; j < columnas; j++) {
            aux = matriz[i][j];
            matriz[i][j] = matriz[j][i];
            matriz[j][i] = aux;
        }
    }
}

//------------------------------------------------------------------------------
// void copiarMatriz6x6(double matriz1[6][6], double matriz2[6][6])
//------------------------------------------------------------------------------
/**
*
* This function copies the values of a 6x6 matrix into another 6x6 matrix.
*
* @param matriz1 The matrix from which to copy the values.
* @param matriz2 The matrix to which the values will be copied.
*/ 
//------------------------------------------------------------------------------

void copiarMatriz6x6(double matriz1[6][6], double matriz2[6][6]){
	for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            matriz2[i][j] = matriz1[i][j];
        }
    }
}

//------------------------------------------------------------------------------
// bool matricesIguales6x6(double matriz1[6][6], double matriz2[6][6], int filas, int columnas)
//------------------------------------------------------------------------------
/**
*
* This function compares the values of two 6x6 matrices and checks if they are equal within a tolerance.
*
* @param matriz1 The first matrix to compare.
* @param matriz2 The second matrix to compare.
* @param filas The number of rows in the matrices.
* @param columnas The number of columns in the matrices.
* @return True if the matrices are equal within the tolerance, false otherwise.
*/ 
//------------------------------------------------------------------------------

bool matricesIguales6x6(double matriz1[6][6], double matriz2[6][6], int filas, int columnas) {
    
    for (int i = 0; i < filas; i++) {
        for (int j = 0; j < columnas; j++) {
            if (fabs(matriz1[i][j] - matriz2[i][j]) > pow(10,-6)) {
                return false;
            }
        }
    }
    return true;
}

//------------------------------------------------------------------------------
// void multMatrices1x6_6x6(double matriz1[1][6], double matriz2[6][6], double result[1][6])
//------------------------------------------------------------------------------
/**
*
* This function performs the multiplication of a 1x6 matrix by a 6x6 matrix and stores the result in a 1x6 matrix.
*
* @param matriz1 The 1x6 matrix.
* @param matriz2 The 6x6 matrix.
* @param result The resulting 1x6 matrix.
*/ 
//------------------------------------------------------------------------------

void multMatrices1x6_6x6(double matriz1[1][6], double matriz2[6][6], double result[1][6]) {
    
	for (int j = 0; j < 6; j++) {
		result[0][j] = 0.0;
		for (int k = 0; k < 6; k++) {
			result[0][j] += matriz1[0][k] * matriz2[k][j];
		}
	}
    
}

//------------------------------------------------------------------------------
// void transpuesta1x6(double matriz[1][6], double result[6][1])
//------------------------------------------------------------------------------
/**
*
* This function transposes a 1x6 matrix by converting it into a 6x1 matrix and storing the result in a 6x1 matrix.
*
* @param matriz The 1x6 matrix to be transposed.
* @param result The resulting 6x1 matrix.
*/ 
//------------------------------------------------------------------------------

void transpuesta1x6(double matriz[1][6], double result[6][1]) {
    for (int i = 0; i < 6; ++i) {
        result[i][0] = matriz[0][i];
    }
}

//------------------------------------------------------------------------------
// double multMatrices1x6_6x1(double matriz1[1][6], double matriz2[6][1])
//------------------------------------------------------------------------------
/**
*
* This function performs matrix multiplication between a 1x6 matrix and a 6x1 matrix. 
* It computes the dot product of the corresponding elements of the two matrices and returns the result.
*
* @param matriz1 The 1x6 matrix.
* @param matriz2 The 6x1 matrix.
* @return The result of the matrix multiplication as a double value.
*/ 
//------------------------------------------------------------------------------

double multMatrices1x6_6x1(double matriz1[1][6], double matriz2[6][1]) {
    double result = 0.0;
    for (int i = 0; i < 6; ++i) {
        result += matriz1[0][i] * matriz2[i][0];
    }
    return result;
}

//------------------------------------------------------------------------------
// void multMatrices6x6_6x1(double matriz1[6][6], double matriz2[6][1], double result[6][1])
//------------------------------------------------------------------------------
/**
*
* This function performs matrix multiplication between a 6x6 matrix and a 6x1 matrix. It computes the dot product of each row of the 6x6 matrix with the corresponding element in the 6x1 matrix and stores the result in a 6x1 matrix.
*
* @param matriz1 The 6x6 matrix.
* @param matriz2 The 6x1 matrix.
* @param result The resulting 6x1 matrix to store the multiplication result.
*/ 
//------------------------------------------------------------------------------

void multMatrices6x6_6x1(double matriz1[6][6], double matriz2[6][1], double result[6][1]) {
 
	for (int i = 0; i < 6; i++) {
		result[i][0] = 0.0;
		for (int k = 0; k < 6; k++) {
			result[i][0] += matriz1[i][k] * matriz2[k][0];
		}
    }
    
}

//------------------------------------------------------------------------------
// void eye6x6(double matriz[6][6])
//------------------------------------------------------------------------------
/**
*
* This function fills a 6x6 matrix with zeros, except for the diagonal elements which are set to 1.0, resulting in an identity matrix.
*
* @param matriz The 6x6 matrix to be generated as an identity matrix.
*/ 
//------------------------------------------------------------------------------

void eye6x6(double matriz[6][6]) {
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            if (i == j) {
                matriz[i][j] = 1.0;
            } else {
                matriz[i][j] = 0.0;
            }
        }
    }
}

//------------------------------------------------------------------------------
// void multMatrices6x1_1x6(double matriz1[6][1], double matriz2[1][6], double result[6][6])
//------------------------------------------------------------------------------
/**
*
* This function performs element-wise multiplication between the corresponding elements of the input matrices, where the first matrix is a 6x1 matrix and the second matrix is a 1x6 matrix. The resulting matrix is a 6x6 matrix.
*
* @param matriz1 The 6x1 matrix to be multiplied.
* @param matriz2 The 1x6 matrix to be multiplied.
* @param result The resulting 6x6 matrix.
*/ 
//------------------------------------------------------------------------------

void multMatrices6x1_1x6(double matriz1[6][1], double matriz2[1][6], double result[6][6]) {
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            result[i][j] = matriz1[i][0] * matriz2[0][j];
        }
    }
}

//------------------------------------------------------------------------------
// void multMatrices6x6(double matriz1[6][6], double matriz2[6][6], double result[6][6])
//------------------------------------------------------------------------------
/**
*
* This function performs matrix multiplication between the given 6x6 matrices, `matriz1` and `matriz2`, and stores the result in the `result` matrix.
*
* @param matriz1 The first 6x6 matrix to be multiplied.
* @param matriz2 The second 6x6 matrix to be multiplied.
* @param result The resulting 6x6 matrix.
*/ 
//------------------------------------------------------------------------------

void multMatrices6x6(double matriz1[6][6], double matriz2[6][6], double result[6][6]) {
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            result[i][j] = 0.0;
            for (int k = 0; k < 6; k++) {
                result[i][j] += matriz1[i][k] * matriz2[k][j];
            }
        }
    }
}

//------------------------------------------------------------------------------
// bool matricesIguales6x1(double matriz1[6][1], double matriz2[6][1])
//------------------------------------------------------------------------------
/**
*
* This function compares each element of the given 6x1 matrices, `matriz1` and `matriz2`, and checks if their absolute difference is within a certain tolerance. If all elements are within the tolerance, the matrices are considered equal and the function returns true; otherwise, it returns false.
*
* @param matriz1 The first 6x1 matrix to be compared.
* @param matriz2 The second 6x1 matrix to be compared.
* @return true if the matrices are equal within the specified tolerance, false otherwise.
*/ 
//------------------------------------------------------------------------------

bool matricesIguales6x1(double matriz1[6][1], double matriz2[6][1]) {
    
    for (int i = 0; i < 6; i++) {
		if (fabs(matriz1[i][0] - matriz2[i][0]) > pow(10,-8)) {
			return false;
		}
    }
    return true;
}

//------------------------------------------------------------------------------
// bool matricesIguales42x1(double matriz1[42][1], double matriz2[42][1])
//------------------------------------------------------------------------------
/**
*
* This function compares each element of the given 42x1 matrices, `matriz1` and `matriz2`, and checks if their absolute difference is within a certain tolerance. If all elements are within the tolerance, the matrices are considered equal and the function returns true; otherwise, it returns false.
*
* @param matriz1 The first 42x1 matrix to be compared.
* @param matriz2 The second 42x1 matrix to be compared.
* @return true if the matrices are equal within the specified tolerance, false otherwise.
*/
//------------------------------------------------------------------------------

bool matricesIguales42x1(double matriz1[42][1], double matriz2[42][1]) {
    
    for (int i = 0; i < 42; i++) {
		if (fabs(matriz1[i][0] - matriz2[i][0]) > pow(10,-9)) {
			return false;
		}
    }
    return true;
}

//------------------------------------------------------------------------------
// void mult1x3(double matriz1[1][3], double matriz2[3][3], double resultado[1][3])
//------------------------------------------------------------------------------
/**
*
* This function performs matrix multiplication between a 1x3 matrix, `matriz1`, and a 3x3 matrix, `matriz2`. The result is a 1x3 matrix, `resultado`, which stores the product of the matrices.
*
* @param matriz1 The 1x3 matrix to be multiplied.
* @param matriz2 The 3x3 matrix to be multiplied.
* @param resultado The resulting 1x3 matrix.
*/
//------------------------------------------------------------------------------

void mult1x3(double matriz1[1][3], double matriz2[3][3], double resultado[1][3]) {
    for (int i = 0; i < 1; i++) {
        for (int j = 0; j < 3; j++) {
			resultado[i][j] = 0;
			for (int k = 0; k < 3; k++) {
				resultado[i][j] += matriz1[i][k] * matriz2[k][j];
			}
		}
        
    }
}

//------------------------------------------------------------------------------
// void vectorAMatriz2(double v[3], double m[1][3])
//------------------------------------------------------------------------------
/**
*
* This function takes a 1D vector `v` of size 3 and converts it into a 1x3 matrix `m`. The elements of the vector are copied into the matrix row-wise.
*
* @param v The 1D vector to be converted.
* @param m The resulting 1x3 matrix.
*/
//------------------------------------------------------------------------------

void vectorAMatriz2(double v[3], double m[1][3]){
	for (int i = 0; i < 3; i++) {
        m[0][i] = v[i];
    }
}

