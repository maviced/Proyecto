#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "include/globales.h"
#include "include/Mjday.h"
#include "include/SAT_Const.h"
#include "include/Position.h"
#include "include/anglesdr.h"
#include "include/IERS.h"
#include "include/timediff.h"
#include "include/DEInteg.h"
#include "include/Accel.h"
#include "include/LTC.h"
#include "include/VarEqn.h"
#include "include/gmst.h"
#include "include/R_z.h"
#include "include/matriz.h"
#include "include/TimeUpdate.h"
#include "include/MeasUpdate.h"
#include "include/AzElPa.h"
#include "include/vector.h"

//------------------------------------------------------------------------------
// int main()
//------------------------------------------------------------------------------
/**
*
* @file EKF_GEOS3.cpp
*
* Initial Orbit Determination using Double-R-Iteration and Extended Kalman Filter methods
*
* Last modified:   2015/08/12   M. Mahooti
*
* References:
*   O. Montenbruck, E. Gill, "Satellite Orbits - Models, Methods, and
*   Applications", Springer Verlag, Heidelberg, 2000.
*   
*   D. Vallado, "Fundamentals of Astrodynamics and Applications", 
*   4th Edition, 2013.
*
*   G. Seeber, "Satellite Geodesy", 2nd Edition, 2003.
*
*
*/ 
//------------------------------------------------------------------------------ 

using namespace std;

double **eopdata;
double pnm[362][362], dpnm[362][362], Cnm[362][362], Snm[362][362];
Param AuxParam;

int main(){
	
	cout << "Programa principal:" << endl;
	cout << setprecision(16) << endl;
	
    FILE *fp;
    int f, c;
	double aux1, aux2, sigma_range, sigma_az, sigma_el, lat, lon, alt, Rs[3], Mjd1, Mjd2, Mjd3, r2[3],v2[3], Y0_apr[6][1], Mjd0, Mjd_UTC, UT1_UTC, TAI_UTC, 
	x_pole, y_pole, UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC, P[6][6], LT[3][3], yPhi[42][1], Phi[6][6], t, Y_old[6][1], t_old, Mjd_TT, Mjd_UT1, theta,
	U[3][3], aux3[3][1], aux4[3][1], s2[3][1], Azim, Elev, dAds[3], dEds[3], s2_v[3], dAds_m[1][3], aux5[1][3], aux6[1][3], dAdY[1][6], K[6][1], Y2[6][1],
	P2[6][6], aux7[3][1], aux8[3][1], dEds_m[1][3], dEdY[1][6], Y3[6][1], dDds[1][3], dDdY[1][6];

	int n_eqn;
	
	for (int i = 0; i < 42; i++) {
		yPhi[i][0] = 0.0;
    }
	
	for (int i = 0; i < 6; i++) {
		for(int j=0; j<6; j++){
			Phi[i][j] = 0.0;
		}
    }
	

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

	AuxParam.Mjd_TT = 0.0;
	AuxParam.n = 0;
	AuxParam.m = 0;
	AuxParam.Mjd_UTC = 0.0;


    fp = fopen("./data/eop19620101.txt", "r");
    if(fp == NULL){
        printf("Fail open eop19620101.txt file\n");
        exit(EXIT_FAILURE);
    }

    eopdata = (double **) malloc(14 * sizeof(double *));
    if(eopdata == NULL){
        printf("eopdata: memory not allocated\n");
        exit(EXIT_FAILURE);
    }
    
    for(int i = 0; i <= 14; i++){
        eopdata[i] = (double *) malloc(21413 * sizeof(double));
        if(eopdata[i] == NULL){
            printf("eopdata[i]: memory not allocated\n");
            exit(EXIT_FAILURE);
        }
    }

    for(int i = 0; i < 21413; i++){
        fscanf(fp, "%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf", &eopdata[0][i], &eopdata[1][i], &eopdata[2][i], &eopdata[3][i], &eopdata[4][i], &eopdata[5][i], &eopdata[6][i], &eopdata[7][i], &eopdata[8][i], &eopdata[9][i], &eopdata[10][i], &eopdata[11][i], &eopdata[12][i]);
    }
    fclose(fp);

    //cout << eopdata[3][2] << endl;
    
	int nobs = 18;
    double obs[nobs][4];

    // read observations
	string tline;
	int Y, M, D, h, m;
    double s, az, el, Dist;
	
	
	ifstream fichGEO("./data/GEOS3.txt");
	
	if (!fichGEO.is_open()) {
		cout << "No se pudo abrir el archivo GEOS3.txt" << endl;
		return 0;
	}
	
	for(int i=0; i<nobs; i++){
		if(!getline(fichGEO, tline)){
			break;
		}
		
		Y = stoi(tline.substr(0, 4));
		M = stoi(tline.substr(5, 2));
		D = stoi(tline.substr(8, 2));
		h = stoi(tline.substr(12, 2));
		m = stoi(tline.substr(15, 2));
		s = stod(tline.substr(18, 6));
		az = stod(tline.substr(25, 8));
		el = stod(tline.substr(35, 7));
		Dist = stod(tline.substr(44, 10));
		obs[i][0] = Mjday(Y,M,D,h,m,s);
		obs[i][1] = Rad*az;
		obs[i][2] = Rad*el;
		obs[i][3] = 1e3*Dist;
		
		
		
		//cout << Y << " " << M << " " << D << " " << h << " " << m << " " << fixed << setprecision(2) << s << " " << fixed << setprecision(4) << az << " " << el << " " << Dist << std::endl;
	}
	
	
	fichGEO.close();
	/*
	for (int i = 0; i < 18; i++) {
        for (int j = 0; j < 4; j++) {
            cout << obs[i][j] << " ";
        }
        cout << endl;
    }*/
	
	sigma_range = 92.5;     // [m]
	sigma_az = 0.0224*Rad;  // [rad]
	sigma_el = 0.0139*Rad;  // [rad]

	// Kaena Point station
	lat = Rad*21.5748;     // [rad]
	lon = Rad*(-158.2706); // [rad]
	alt = 300.20;          // [m]
	
	Position(lon, lat, alt, Rs);
	/*
	cout << Rs[0] << endl;
	cout << Rs[1] << endl;
	cout << Rs[2] << endl;*/
	
	Mjd1 = obs[0][0];
	Mjd2 = obs[8][0];
	Mjd3 = obs[17][0];
	/*
	cout << Mjd1 << endl;
	cout << Mjd2 << endl;
	cout << Mjd3 << endl;*/
	
	anglesdr ( obs[0][1],obs[8][1],obs[17][1],obs[0][2],obs[8][2],obs[17][2],Mjd1,Mjd2,Mjd3,Rs,Rs,Rs, r2,v2 );
	
	//Y0_apr = [r2;v2]
	Y0_apr[0][0] = r2[0];
	Y0_apr[1][0] = r2[1];
	Y0_apr[2][0] = r2[2];
	Y0_apr[3][0] = v2[0];
	Y0_apr[4][0] = v2[1];
	Y0_apr[5][0] = v2[2];
	
	/*
	for(int i=0; i<6; i++){
		cout << Y0_apr[i][0] << endl;
	}*/
	
	Mjd0 = Mjday(1995,1,29, 02,38,00.0);
	
	Mjd_UTC = obs[8][0];
	IERS (eopdata, Mjd_UTC, &UT1_UTC, &TAI_UTC, &x_pole, &y_pole);
	timediff(UT1_UTC, TAI_UTC, &UT1_TAI, &UTC_GPS, &UT1_GPS, &TT_UTC, &GPS_UTC);
	
	/*
	cout << Mjd0 << endl;
	cout << Mjd_UTC << endl;
	cout << UT1_UTC << endl;
	cout <<  TAI_UTC << endl;
	cout <<  x_pole << endl;
	cout <<  y_pole << endl;
	cout << UT1_TAI  << endl;
	cout << UTC_GPS  << endl;
	cout <<  UT1_GPS << endl;
	cout <<  TT_UTC << endl;
	cout <<  GPS_UTC << endl;*/
	
	AuxParam.Mjd_TT  = Mjd_UTC + TT_UTC/86400;
	AuxParam.n       = 10;
	AuxParam.m       = 10;

	n_eqn  = 6;
	
	DEInteg (&Accel,0,-(obs[8][0]-Mjd0)*86400.0,1e-13,1e-6,6,Y0_apr);
	/*
	for (int i = 0; i < 6; i++) {
        cout << Y0_apr[i][0] << endl;
    }*/
	
	for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            P[i][j] = 0.0;
        }
    }
  
  
	for (int i=0; i<3; i++){
		P[i][i]=1e8;
	}
	for (int i=3; i<6; i++){
		P[i][i]=1e3;
	}
	
	LTC(lon,lat, LT);
	/*
	for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            cout << LT[i][j] << " ";
        }
        cout << endl;
    }*/
	

	// Measurement loop
	t = 0.0;
	
	for (int i=0; i<nobs; i++){    
		// Previous step
		t_old = t;
		for(int j=0; j<6; j++){
			Y_old[j][0] = Y0_apr[j][0];
		}
		
		// Time increment and propagation
		Mjd_UTC = obs[i][0];                      // Modified Julian Date
		t       = (Mjd_UTC-Mjd0)*86400.0;         // Time since epoch [s]
		
		IERS (eopdata, Mjd_UTC, &UT1_UTC, &TAI_UTC, &x_pole, &y_pole);
		timediff(UT1_UTC, TAI_UTC, &UT1_TAI, &UTC_GPS, &UT1_GPS, &TT_UTC, &GPS_UTC);
		Mjd_TT = Mjd_UTC + TT_UTC/86400;
		Mjd_UT1 = Mjd_TT + (UT1_UTC-TT_UTC)/86400.0;
		AuxParam.Mjd_UTC = Mjd_UTC;
		AuxParam.Mjd_TT = Mjd_TT;
			
		for (int ii=1; ii<7; ii++){
			yPhi[ii-1][0] = Y_old[ii-1][0];
			for (int j=1; j<7; j++){  
				if (ii==j){ 
					yPhi[6*j+ii-1][0] = 1; 
				}else{
					yPhi[6*j+ii-1][0] = 0;
				}
			}
		}
		
		
		DEInteg (&VarEqn,0,t-t_old,1e-13,1e-6,42,yPhi);
		
		
		// Extract state transition matrices
		for (int j=1; j<7; j++){
			for(int k=0; k<6; k++){
				Phi[k][j-1] = yPhi[6*j+k][0];
			}
		}
		/*
		for (int i = 0; i < 6; i++) {
			for (int j = 0; j < 6; j++) {
				cout << Phi[i][j] << " ";
			}
			cout << endl;
		}*/
		
		DEInteg (&Accel,0,t-t_old,1e-13,1e-6,6,Y_old);
		
		// Topocentric coordinates
		theta = gmst(Mjd_UT1);                    // Earth rotation
		R_z(theta, U);
		double r[3][1] = {{Y_old[0][0]}, {Y_old[1][0]}, {Y_old[2][0]}};
		mult3x1(U,r, aux3);
		for(int j=0; j<3; j++){
			aux4[j][0] = aux3[j][0] - Rs[j];
		}
		mult3x1(LT, aux4, s2);					  // Topocentric position [m]
		/*
		for(int j=0; j<3; j++){
			cout << s2[j][0] << endl;
		}*/
		
		// Time update
		TimeUpdate(P, Phi);
		/*
		for (int i = 0; i < 6; i++) {
			for (int j = 0; j < 6; j++) {
				cout << P[i][j] << " ";
			}
			cout << endl;
		}*/
			
		// Azimuth and partials
		matrizAVector(s2, s2_v);
		AzElPa(s2_v, &Azim, &Elev, dAds, dEds);     // Azimuth, Elevation
		
		vectorAMatriz2(dAds, dAds_m);
		mult1x3(dAds_m,LT, aux5);
		mult1x3(aux5, U, aux6);
		/*
		for(int j=0; j<3; j++){
			cout <<aux6[0][j] << endl;
		}*/
		
		dAdY[0][0] = aux6[0][0];
		dAdY[0][1] = aux6[0][1];
		dAdY[0][2] = aux6[0][2];
		dAdY[0][3] = 0.0;
		dAdY[0][4] = 0.0;
		dAdY[0][5] = 0.0;
		
		// Measurement update
		MeasUpdate ( Y_old, obs[i][1], Azim, sigma_az, dAdY, P, 6, K, Y2, P2 );
		/*for(int j=0; j<6; j++){
			cout << Y2[j][0] << endl;
		}
		for (int i = 0; i < 6; i++) {
			for (int j = 0; j < 6; j++) {
				cout << P2[i][j] << " ";
			}
			cout << endl;
		}*/
		
		
		// Elevation and partials
		r[0][0] = Y2[0][0];
		r[1][0] = Y2[1][0];
		r[2][0] = Y2[2][0];
		/*
		for(int j=0; j<3; j++){
			cout << r[j][0] << endl;
		}*/
		
		mult3x1(U, r, aux7);
		for(int j=0; j<3; j++){
			aux8[j][0] = aux7[j][0] - Rs[j];
		}
		mult3x1(LT, aux8, s2);					  // Topocentric position [m]
		/*for(int j=0; j<3; j++){
			cout << s2[j][0] << endl;
		}*/
		
		matrizAVector(s2, s2_v);
		AzElPa(s2_v, &Azim, &Elev, dAds, dEds);     // Azimuth, Elevation
		
		
		vectorAMatriz2(dEds, dEds_m);
		mult1x3(dEds_m,LT, aux5);
		mult1x3(aux5, U, aux6);
		/*
		for(int j=0; j<3; j++){
			cout <<aux6[0][j] << endl;
		}*/
		
		dEdY[0][0] = aux6[0][0];
		dEdY[0][1] = aux6[0][1];
		dEdY[0][2] = aux6[0][2];
		dEdY[0][3] = 0.0;
		dEdY[0][4] = 0.0;
		dEdY[0][5] = 0.0;
		
		
		// Measurement update
		MeasUpdate ( Y2, obs[i][2], Elev, sigma_el, dEdY, P2, 6, K, Y3, P );
		/*for(int j=0; j<6; j++){
			cout << Y3[j][0] << endl;
		}
		for (int i = 0; i < 6; i++) {
			for (int j = 0; j < 6; j++) {
				cout << P[i][j] << " ";
			}
			cout << endl;
		}*/

		
		// Range and partials
		r[0][0] = Y3[0][0];
		r[1][0] = Y3[1][0];
		r[2][0] = Y3[2][0];
		
		mult3x1(U, r, aux7);
		for(int j=0; j<3; j++){
			aux8[j][0] = aux7[j][0] - Rs[j];
		}
		mult3x1(LT, aux8, s2);					  // Topocentric position [m]
		/*for(int j=0; j<3; j++){
			cout << s2[j][0] << endl;
		}*/
		
		matrizAVector(s2, s2_v);
		Dist = norm(s2_v, 3);
		//cout << Dist << endl;
		
		for(int j=0; j<3; j++){
			dDds[0][j] = s2_v[j]/Dist;
		}
		/*for(int j=0; j<3; j++){
			cout << dDds[0][j] << endl;
		}*/
		
		mult1x3(dDds,LT, aux5);
		mult1x3(aux5, U, aux6);
		/*
		for(int j=0; j<3; j++){
			cout <<aux6[0][j] << endl;
		}*/
		
		dDdY[0][0] = aux6[0][0];
		dDdY[0][1] = aux6[0][1];
		dDdY[0][2] = aux6[0][2];
		dDdY[0][3] = 0.0;
		dDdY[0][4] = 0.0;
		dDdY[0][5] = 0.0;
		
		/*for(int j=0; j<6; j++){
			cout << dDdY[0][j] << endl;
		}*/
		
		// Measurement update
		MeasUpdate ( Y3, obs[i][3], Dist, sigma_range, dDdY, P, 6 , K, Y0_apr, P2);
		for (int i = 0; i < 6; i++) {
			for (int j = 0; j < 6; j++) {
				P[i][j] = P2[i][j];
			}
		}
		/*for(int j=0; j<6; j++){
			cout << Y0_apr[j][0] << endl;
		}
		for (int i = 0; i < 6; i++) {
			for (int j = 0; j < 6; j++) {
				cout << P[i][j] << " ";
			}
			cout << endl;
		}*/
	}
	
	
	IERS (eopdata, obs[17][0], &UT1_UTC, &TAI_UTC, &x_pole, &y_pole);
	timediff(UT1_UTC, TAI_UTC, &UT1_TAI, &UTC_GPS, &UT1_GPS, &TT_UTC, &GPS_UTC);
	Mjd_TT = Mjd_UTC + TT_UTC/86400;
	AuxParam.Mjd_UTC = Mjd_UTC;
	AuxParam.Mjd_TT = Mjd_TT;
	
	//cout << AuxParam.Mjd_UTC << endl;
	//cout << AuxParam.Mjd_TT << endl;
	
	DEInteg (&Accel,0,-(obs[17][0]-obs[0][0])*86400.0,1e-13,1e-6,6,Y0_apr);

	/*for(int j=0; j<6; j++){
		cout << Y0_apr[j][0] << endl;
	}*/


	double Y_true[6][1] = {{5753.173e3}, {2673.361e3}, {3440.304e3}, {4.324207e3}, {-1.924299e3}, {-5.728216e3}};
	

	printf("\nError of Position Estimation\n");
	printf("dX%10.1f [m]\n",Y0_apr[0][0]-Y_true[0][0]);
	printf("dY%10.1f [m]\n",Y0_apr[1][0]-Y_true[1][0]);
	printf("dZ%10.1f [m]\n",Y0_apr[2][0]-Y_true[2][0]);
	printf("\nError of Velocity Estimation\n");
	printf("dVx%8.1f [m/s]\n",Y0_apr[3][0]-Y_true[3][0]);
	printf("dVy%8.1f [m/s]\n",Y0_apr[4][0]-Y_true[4][0]);
	printf("dVz%8.1f [m/s]\n",Y0_apr[5][0]-Y_true[5][0]);
	

    
    for(int i = 0; i <= 14; i++){
        free(eopdata[i]);
    }
    free(eopdata);
    
    return 0;

}
	