#pragma once
#include "alltypes.h"
#include <vector>
#include <array>
#include <fstream>
#include <sstream>
#include <string>
#include <iostream>






    
    double* Cnm = nullptr;
    double* Snm = nullptr;
    int** Order2 = nullptr;

   double PfBel(double n, double m) {
        if (m == 0) return std::sqrt(2. * n + 1);
        double e = (2. * n + 1.) * 2.;
        double c = n;
        double a = (n + m) / 4.;
        double ac = a / c;
        for (int i = 1; i <= (int)m - 1; i++) {
            ac *= (n + m - i) / 4. / (n - i);
        }
        return std::sqrt(e * ac);
    }

   int order2(int n, int m) {
        if (n == 0) return m;
        return order2(n - 1, n - 1) + m + 1;
    }

  void setOrder2(int nmax) {
        Order2 = new int* [nmax + 1];
        for (int i = 0; i <= nmax; i++) {
            Order2[i] = new int[nmax + 1];
            for (int j = 0; j <= i; j++)
                Order2[i][j] = order2(i, j);
        }
    }

  void importStokesBelikov(const std::string& path, int nmax) {
        unsigned int Coff_Dim = (nmax + 2) * (nmax + 1) / 2;
        Cnm = (double*)calloc(Coff_Dim, sizeof(double));
        Snm = (double*)calloc(Coff_Dim, sizeof(double));
        if (!Cnm || !Snm) {
            std::cerr << "fatal: out of memory (fillStokes).\n";
            exit(EXIT_FAILURE);
        }

        setOrder2(nmax);

        std::ifstream file(path);
        if (!file) {
            std::cerr << "Error: can't open Stokes file " << path << "\n";
            exit(EXIT_FAILURE);
        }

        std::string line;
        int count = 0;
        while (std::getline(file, line)) {
            if (line.empty() || line[0] == '#') continue;
            std::istringstream iss(line);
    
            int n, m;
            double C, S;
           
            iss >> n >> m >> C >> S;
            if (n > nmax || m > n) continue;

            double norm = PfBel(n, m);
            int index = order2(n, m);
            Cnm[index] = C * norm;
            Snm[index] = S * norm;

        }

        setOrder2(nmax);
    }

 
   void gravityBelikov(double r, double lat, double lon, int nmax, std::array<double, 3>& Result)
  {
       
       double lat_rad = lat * M_PI / 180.0;  
       double lon_rad = lon * M_PI / 180.0;   


     double R = EARTH_RADIUS / r;
     double theta = M_PI / 2.0 - lat_rad; 
     double t = std::cos(theta);            
     double u = std::sin(theta);                  

        
        int size = nmax + 1;
        double mu = EARTH_MU / r;

        double Vr = 0., Vt = 0., Vtl = 0.;
        double** P = new double* [size];
        for (int i = 0; i < size; i++) P[i] = new double[size + 3]();

        P[0][0] = 1.;
        
        double cosl = std::cos(lon_rad);       
        double sinl = std::sin(lon_rad);      
        P[0][size + 1] = 1.;
        P[0][size + 2] = 0.;

        for (int n = 1; n < size; n++) {
            for (int m = 0; m <= n; m++) {
                double a = (n - 1 >= m) ? P[n - 1][m] : 0.;
                double b = (n - 1 >= m + 1) ? P[n - 1][m + 1] : 0.;
                double c = (m == 0 && n - 1 >= 1) ? -0.25 * P[n - 1][1] :
                    (m > 0 && n - 1 >= m - 1) ? P[n - 1][m - 1] : 0.;

                P[n][m] = t * a - u / 4. * (b - 4. * c);
                P[m][n + 1] = -n * (u * a + t / 4. * (b - 4. * c));
            }
            P[n][size + 1] = cosl * P[n - 1][size + 1] - sinl * P[n - 1][size + 2];
            P[n][size + 2] = cosl * P[n - 1][size + 2] + sinl * P[n - 1][size + 1];
        }

        for (int m = 0; m < size; m++) {
            double cl = P[m][size + 1];
            double sl = P[m][size + 2];
            for (int n = m; n < size; n++) {
                int ind = Order2[n][m];
                double cnm = Cnm[ind];
                double snm = Snm[ind];
                double tt = std::pow(R, n) * P[n][m];

                Vr += mu * (n + 1) / r * (-cnm * tt * cl - snm * tt * sl);
                Vt += mu * std::pow(R, n) * P[m][n + 1] * (cnm * cl + snm * sl);
                Vtl += mu * m * (tt / u) * (snm * cl - cnm * sl);
            }
        }



       
        
        
       
        double ax = Vr * u * cosl + Vt * t * cosl / r - Vtl * sinl / r;
        double ay = Vr * u * sinl + Vt * t * sinl / r + Vtl * cosl / r;
        double az = Vr * t - Vt * u / r;



       

        Result[0] = ax;
        Result[1] = ay;
        Result[2] = az;


        for (int i = 0; i < size; i++) delete[] P[i];
        delete[] P;

      
    }

