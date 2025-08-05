#pragma once
#include <vector>
#include <cmath>
#include <iostream>
#define _USE_MATH_DEFINES
#include <map> 
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cstring>



#define EARTH_RADIUS 6378136.30 
#define EARTH_MU 398600441500000.0 
	
	void gravityBelikov(double r, double lat, double lon, int n, std::array<double, 3>& Result);
	void gravityCunningham(double r, double lat, double lon, int n, std::array<double, 3>& Result);


std::map<int, std::string> gravityModels = {
	{1, "EGM96.dat"},
	{2, "egm2008.dat"}
};

int selectedModel = 1; 

void selectGravityModel(int& selectedModel, const std::map<int, std::string>& models) {
	std::cout << "Available gravity models:\n";
	for (const auto& pair : models) {
		std::cout << " " << pair.first << ". " << pair.second << "\n";
	}
	std::cout << "Enter model number: ";
	int choice;
	std::cin >> choice;
	if (models.find(choice) != models.end()) {
		selectedModel = choice;
		std::cout << "Selected model: " << models.at(choice) << "\n";
	}
	else {
		std::cout << "Invalid model number.\n";
	}
}





double* Cnm = nullptr;
double* Snm = nullptr;
double* Cnm1 = nullptr;
double* Snm1 = nullptr;
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






void importStokesCombined(const std::string& path, int nmax) {
  
    unsigned int Coff_Dim = (nmax + 2) * (nmax + 1) / 2;


    Cnm = (double*)calloc(Coff_Dim, sizeof(double));
    Snm = (double*)calloc(Coff_Dim, sizeof(double));
    Cnm1 = (double*)calloc(Coff_Dim, sizeof(double));
    Snm1 = (double*)calloc(Coff_Dim, sizeof(double));

    if (!Cnm || !Snm || !Cnm1 || !Snm1) {
        std::cerr << "fatal: out of memory (importStokesCombined).\n";
        exit(EXIT_FAILURE);
    }

    setOrder2(nmax);

    std::ifstream file(path);
    if (!file) {
        std::cerr << "Error: can't open Stokes file " << path << "\n";
        exit(EXIT_FAILURE);
    }

    std::string line;
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;

        std::istringstream iss(line);
        int n, m;
        double C, S;

   
        iss >> n >> m >> C >> S;
        if (iss.fail()) continue;

        if (n > nmax || m > n) continue;

        int idx = order2(n, m);
        double norm = PfBel(n, m);

      
        Cnm[idx] = C * norm;
        Snm[idx] = S * norm;

       
        Cnm1[idx] = C;
        Snm1[idx] = S;
    }
}



void freeStokes(int nmax) {
	if (Cnm)   free(Cnm);
	if (Snm)   free(Snm);
	if (Cnm1)  free(Cnm1);
	if (Snm1)  free(Snm1);

	Cnm = Snm = Cnm1 = Snm1 = nullptr;

	if (Order2) {
		for (int i = 0; i <= nmax; i++) {
			if (Order2[i]) free(Order2[i]);
		}
		free(Order2);
		Order2 = nullptr;
	}
}
