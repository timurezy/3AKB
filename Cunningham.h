
#pragma once

#include <vector>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include "alltypes.h"



	void ManualLegendre(double P[], double Pn_1[], double Pn_2[], double sinphi, double cosphi) {


		P[0] = 1;									// Algorithm Seed
		Pn_1[0] = P[0];

		P[0] = sinphi * Pn_1[0] * sqrt(3.0);		// Legendre Polynomial P(n,m) = P(1,0);
		P[1] = cosphi * Pn_1[0] * sqrt(3.0);		// Legendre Polynomial P(n,m) = P(1,1);

		// Copy Back Results

		Pn_2[0] = Pn_1[0];

		Pn_1[0] = P[0];
		Pn_1[1] = P[1];


	}


	void CopyBackRecursion(double P[], double Pn_1[], double Pn_2[], int n) {

		int i = 0;

		for (i = 0; i <= n; i++) {

			Pn_2[i] = Pn_1[i];
			Pn_1[i] = P[i];

		}



	}


	inline void Normalized_Legendre(double P[], double Pn_1[], double Pn_2[], double sinphi, double cosphi, int n) {

		//   Recursive Computation Of The Associated Legendre Polynomials of Order n

		//   Author: Diego Garcнa Pardo (UNIVERSITY CARLOS III OF MADRID)
		//   
		//   
		//   Requires as Input lower order Legendre Polynomials P @ n-1 and P @ n-2
		//   
		//   Also requires the introduction of the values of the sine and cosine of
		//   the longitude angle (measured from the x-y plane (horizontal with
		//   z-axis over the north pole)
		//   
		//   Therefore this function is ONLY VALID for greater or equal than 2 polynomial

		//  ONLY VALID FOR n > 2

		int i = 0;
		int* m;
		double chi_1, chi_2, chi_3, chi_4;

		m = (int*)calloc((n + 1), sizeof(int));
		if (!m) {
			
			std::cout << "fatal: out of memory (m).\n";
			exit(EXIT_FAILURE);
		}

		// Create Vector m (ALF)
		for (i = 0; i <= n; i++) {
			m[i] = i;
		}


		for (i = 0; i <= (n - 2); i++) {

			// Normalizing Parameters
			chi_1 = sqrt((double)(2 * n + 1) * (n - m[i]) / ((2 * n - 1) * (n + m[i])));
			chi_2 = sqrt((double)(2 * n + 1) * (n - m[i]) * (n - m[i] - 1) / ((2 * n - 3) * (n + m[i]) * (n + m[i] - 1)));

			// Legendre Polynomial Value
			P[i] = (double)1 / (n - m[i]) * ((2 * n - 1) * sinphi * Pn_1[i] * chi_1 - (n + m[i] - 1) * Pn_2[i] * chi_2);

		}

		// Additional Normalizing Parameters
		chi_3 = sqrt((double)2 * n + 1);
		chi_4 = sqrt((double)(2 * n + 1) / (2 * n));

		// Resting Legendre Polynomial Values
		P[n - 1] = (double)sinphi * Pn_1[n - 1] * chi_3;
		P[n] = (double)cosphi * Pn_1[n - 1] * chi_4;

		CopyBackRecursion(P, Pn_1, Pn_2, n);

		free(m);
		m = NULL;
	}


	void Acc_Field(double GM, double R, double r, double phi, double lambda, double* Cnm, double* Snm, unsigned int Model_Size, double* _output) {


		// Function Local Variables;
		//double* rdotdot;
		double* V, * W;
		double NF1, NF2, NF3;
		double* CC, * SS;
		double d2x = 0;
		double d2y = 0;
		double d2z = 0;
		unsigned int n, N, m, M;
		unsigned int i = 0;
		unsigned  int I = 0;
		double* P, * Pn_1, * Pn_2;

		double sinphi, cosphi;


		V = (double*)calloc(Model_Size, sizeof(double));
		W = (double*)calloc(Model_Size, sizeof(double));

		CC = (double*)calloc(Model_Size, sizeof(double));
		SS = (double*)calloc(Model_Size, sizeof(double));

		sinphi = sin(phi);
		cosphi = cos(phi);


		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		// Allocating Legendre Vectors



		P = (double*)calloc(Model_Size, sizeof(double));
		Pn_1 = (double*)calloc(Model_Size, sizeof(double));
		Pn_2 = (double*)calloc(Model_Size, sizeof(double));
		if (!V || !W || !CC || !SS || !P || !Pn_1 || !Pn_2) {
		
			std::cout << "fatal: out of memory (Acc_Field).\n";
			exit(EXIT_FAILURE);
		}

		ManualLegendre(P, Pn_1, Pn_2, sinphi, cosphi);
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		// Zero Order Model

		n = 0;									// Model Order Counter
		m = 0;									// Associated Index (Degree)
		N = n + 1;								// Counters for functions V & W (which are 1 order higher)

		for (i = 0; i <= n; i++) {

			CC[i] = Cnm[i + I];
			SS[i] = Snm[i + I];

		}
		I = I + n + 1;							// Vector Counter in Snm And Cnm
		// Zonal Coefficients Associated Acceleration @(m = 0)

		for (M = 0; M <= N; M++) {

			V[M] = pow(R / r, (N + 1)) * Pn_1[M] * cos(M * lambda);
			W[M] = pow(R / r, (N + 1)) * Pn_1[M] * sin(M * lambda);

		}

		NF1 = sqrt((double)0.5 * (2 * n + 1) * (n + 2) * (n + 1) / (2 * n + 3));
		NF2 = sqrt((double)(2 * n + 1) * (n + m + 1) / (2 * n + 3) / (n - m + 1));

		d2x = -CC[m] * V[m + 1] * NF1;
		d2y = -CC[m] * W[m + 1] * NF1;
		d2z = (n + 1) * (-CC[m] * V[m] - SS[m] * W[m]) * NF2;
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		// Counter increase to order 1

		n = n + 1;            // n = 1;
		N = n + 1;

		for (i = 0; i <= n; i++) {

			CC[i] = Cnm[i + I];
			SS[i] = Snm[i + I];

		}
		I = I + n + 1;

		Normalized_Legendre(P, Pn_1, Pn_2, sinphi, cosphi, N);
		for (M = 0; M <= N; M++) {

			V[M] = pow(R / r, (N + 1)) * Pn_1[M] * cos(M * lambda);
			W[M] = pow(R / r, (N + 1)) * Pn_1[M] * sin(M * lambda);

		}

		// Zonal Coefficients Associated Acceleration @(m = 0)
		m = 0;
		NF1 = sqrt((double)0.5 * (2 * n + 1) * (n + 2) * (n + 1) / (2 * n + 3));
		NF2 = sqrt((double)(2 * n + 1) * (n + m + 1) / (2 * n + 3) / (n - m + 1));

		d2x = -CC[m] * V[m + 1] * NF1 + d2x;
		d2y = -CC[m] * W[m + 1] * NF1 + d2y;
		d2z = (n - m + 1) * (-CC[m] * V[m] - SS[m] * W[m]) * NF2 + d2z;

		// Acceleration @(m = 1), Scheme Exception Due To Normalization

		m = 1;
		NF1 = sqrt((double)(2 * n + 1) / (2 * n + 3) * (n + m + 2) * (n + m + 1));
		NF2 = sqrt((double)2 * (2 * n + 1) / (2 * n + 3) / (n - m + 2) / (n - m + 1));
		NF3 = sqrt((double)(2 * n + 1) * (n + m + 1) / (2 * n + 3) / (n - m + 1));



		d2x = (double)0.5 * ((-CC[m] * V[m + 1] - SS[m] * W[m + 1]) * NF1 + (n - m + 2) * (n - m + 1) * (+CC[m] * V[m - 1] + SS[m] * W[m - 1]) * NF2) + d2x;
		d2y = (double)0.5 * ((-CC[m] * W[m + 1] + SS[m] * V[m + 1]) * NF1 + (n - m + 2) * (n - m + 1) * (-CC[m] * V[m - 1] + SS[m] * V[m - 1]) * NF2) + d2y;
		d2z = (double)(n - m + 1) * (-CC[m] * V[m] - SS[m] * W[m]) * NF3 + d2z;

		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		// End Of Scheme Exceptions, LOOP
		//klmClock(true, 0);
		for (n = 2; n < (Model_Size - 1); n++) {

			N = n + 1;


			for (i = 0; i <= n; i++) {

				CC[i] = Cnm[i + I];
				SS[i] = Snm[i + I];

			}

			I = I + n + 1;
			Normalized_Legendre(P, Pn_1, Pn_2, sinphi, cosphi, N);

			for (M = 0; M <= N; M++) {

				V[M] = pow(R / r, (N + 1)) * Pn_1[M] * cos(M * lambda);
				W[M] = pow(R / r, (N + 1)) * Pn_1[M] * sin(M * lambda);

			}

			// Zonal Coefficients Associated Acceleration @(m = 0)
			m = 0;
			NF1 = sqrt((double)0.5 * (2 * n + 1) * (n + 2) * (n + 1) / (2 * n + 3));
			NF2 = sqrt((double)(2 * n + 1) * (n + m + 1) / (2 * n + 3) / (n - m + 1));

			d2x = -CC[m] * V[m + 1] * NF1 + d2x;
			d2y = -CC[m] * W[m + 1] * NF1 + d2y;
			d2z = (n - m + 1) * (-CC[m] * V[m] - SS[m] * W[m]) * NF2 + d2z;

			// Acceleration @(m = 1), Scheme Exception Due To Normalization

			m = 1;
			NF1 = sqrt((double)(2 * n + 1) / (2 * n + 3) * (n + m + 2) * (n + m + 1));
			NF2 = sqrt((double)2 * (2 * n + 1) / (2 * n + 3) / (n - m + 2) / (n - m + 1));
			NF3 = sqrt((double)(2 * n + 1) * (n + m + 1) / (2 * n + 3) / (n - m + 1));

			d2x = (double)0.5 * ((-CC[m] * V[m + 1] - SS[m] * W[m + 1]) * NF1 + (n - m + 2) * (n - m + 1) * (+CC[m] * V[m - 1] + SS[m] * W[m - 1]) * NF2) + d2x;
			d2y = (double)0.5 * ((-CC[m] * W[m + 1] + SS[m] * V[m + 1]) * NF1 + (n - m + 2) * (n - m + 1) * (-CC[m] * V[m - 1] + SS[m] * V[m - 1]) * NF2) + d2y;
			d2z = (double)(n - m + 1) * (-CC[m] * V[m] - SS[m] * W[m]) * NF3 + d2z;


			// Generality of Terms
			for (m = 2; m <= n; m++) {

				NF1 = sqrt((double)(2 * n + 1) / (2 * n + 3) * (n + m + 2) * (n + m + 1));
				NF2 = sqrt((double)(2 * n + 1) / (2 * n + 3) / (n - m + 2) / (n - m + 1));
				NF3 = sqrt((double)(2 * n + 1) * (n + m + 1) / (2 * n + 3) / (n - m + 1));

				d2x = d2x + 0.5 * ((-CC[m] * V[m + 1] - SS[m] * W[m + 1]) * NF1 + (n - m + 2) * (n - m + 1) * (+CC[m] * V[m - 1] + SS[m] * W[m - 1]) * NF2);
				d2y = d2y + 0.5 * ((-CC[m] * W[m + 1] + SS[m] * V[m + 1]) * NF1 + (n - m + 2) * (n - m + 1) * (-CC[m] * V[m - 1] + SS[m] * V[m - 1]) * NF2);
				d2z = d2z + (n - m + 1) * (-CC[m] * V[m] - SS[m] * W[m]) * NF3;


			}
		}

		_output[0] = d2x * GM / pow(R, 2);
		_output[1] = d2y * GM / pow(R, 2);
		_output[2] = d2z * GM / pow(R, 2);

		free(V);
		free(W);

		free(CC);
		free(SS);
		free(P);
		free(Pn_1);
		free(Pn_2);
		V = NULL; W = NULL; CC = NULL; SS = NULL; P = NULL; Pn_1 = NULL; Pn_2 = NULL;
	}







	/*
	inline void importStokesCunningham(const std::string& filename, int nmax) {
		std::ifstream file(filename);
		if (!file) {
			std::cout << "Error opening EGM96.dat file.";
			return;
		}

		int n, m;
		double c, s;
		while (file >> n >> m >> c >> s) {
			int index = n * (n + 1) / 2 + m;
			Cnm1[index] = c;
			Snm1[index] = s;
		}
	}
	*/


	void gravityCunningham(double r, double lat, double lon, int n, std::array<double, 3>& Result)
	{
		
		 
		double phi = lat * M_PI / 180.0;    
		double lambda = lon * M_PI / 180.0;  


		double* rdotdot;


		rdotdot = (double*)calloc(6, sizeof(double));
		if (!rdotdot) {
			std::cout << "fatal: out of memory (rdotdot).\n";
			exit(EXIT_FAILURE);
		};

		

	 Acc_Field(EARTH_MU, EARTH_RADIUS, r, phi, lambda, Cnm1, Snm1, n + 2, rdotdot); 

	
		

		Result[0] = rdotdot[0];
		Result[1] = rdotdot[1];
		Result[2] = rdotdot[2];
		

		free(rdotdot);
		rdotdot = NULL;

	
	}


    