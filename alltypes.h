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
#include <array>
#include <random>
#define M_PI 3.14159265358979323846


#define EARTH_RADIUS 6378136.30 
#define EARTH_MU 398600441500000.0 
	
	void gravityBelikov(double r, double lat, double lon, int n, std::array<double, 3>& Result);
	void gravityCunningham(double r, double lat, double lon, int n, std::array<double, 3>& Result);


std::map<int, std::string> gravityModels = {
	{1, "../EGM96.dat"},
	{2, "../EGM08.dat"}
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
std::vector<std::vector<double>> _c;
std::vector<std::vector<double>> _s;
int importflag = 0;



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

    
    _c = std::vector<std::vector<double>>(nmax + 1);
    _s = std::vector<std::vector<double>>(nmax + 1);
    for (int n = 0; n <= nmax; ++n) {
        _c[n] = std::vector<double>(n + 1, 0.0);
        _s[n] = std::vector<double>(n + 1, 0.0);
    }

    
    std::string line;
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;

        std::istringstream iss(line);
        int n, m;
        double C, S;
        iss >> n >> m >> C >> S;
        if (iss.fail()) continue;
        if (n > nmax || m > n) break;

        int idx = order2(n, m);
        double norm = PfBel(n, m);

        // Çàïèñü â Belikov Cunningham
        Cnm[idx] = C * norm;
        Snm[idx] = S * norm;
        Cnm1[idx] = C;
        Snm1[idx] = S;
     
        // Çàïèñü â _c/_s (gravity_stokes)
        if (m <= n) {
            _c[n][m] = C;
            _s[n][m] = S;
        }
    }
    importflag = 1;
    std::cout << "INFO  | Gravity field loaded successfully into both formats.\n";
}



void freeStokes(int nmax) {
  
    if (Cnm) { free(Cnm);  Cnm = nullptr; }
    if (Snm) { free(Snm);  Snm = nullptr; }
    if (Cnm1) { free(Cnm1); Cnm1 = nullptr; }
    if (Snm1) { free(Snm1); Snm1 = nullptr; }

  
    _c.clear();
    _s.clear();
    _c.shrink_to_fit();
    _s.shrink_to_fit();

 
    if (Order2) {
        for (int i = 0; i <= nmax; ++i) {
            if (Order2[i]) {
                free(Order2[i]);
                Order2[i] = nullptr;
            }
        }
        free(Order2);
        Order2 = nullptr;
    }

    importflag = 0;
}






std::array<double, 3> RLatLonToXYZ(double radius, double latitude, double longitude) {
   
    double latRad = latitude * M_PI / 180.0;
    double lonRad = longitude * M_PI / 180.0;

    double x = radius * std::cos(latRad) * std::cos(lonRad);
    double y = radius * std::cos(latRad) * std::sin(lonRad);
    double z = radius * std::sin(latRad);

    
    return { x, y, z };
}






double first_speed(double radius, double GM = 3.986004418e14) {
    return std::sqrt(GM / radius);
}


std::array<double, 3> random_velocity(double radius) {
    double speed = first_speed(radius);


    static std::random_device rd;
    static std::mt19937 gen(rd());
    static std::uniform_real_distribution<double> dist(0.0, 1.0);

  
    double theta = acos(1.0 - 2.0 * dist(gen)); 
    double phi = 2.0 * M_PI * dist(gen);       

    std::array<double, 3> v;
    v[0] = speed * sin(theta) * cos(phi); 
    v[1] = speed * sin(theta) * sin(phi); 
    v[2] = speed * cos(theta);            

    return v;
}

std::array<double, 3> orbit_velocity(double radius, double latitude, double longitude) {
    double speed = first_speed(radius);

    std::array<double, 3> v;
    v[0] = -1.0 * speed * sin(longitude * M_PI / 180.0);
    v[1] = speed * cos(longitude * M_PI / 180.0);
    v[2] = 0.0;

    return v;
}


std::array<double, 3> XYZtoRLatLon(double x, double y, double z) {
    double radius = std::sqrt(x * x + y * y + z * z);       
    double latitude = std::asin(z / radius) * 180.0 / M_PI;  
    double longitude = std::atan2(y, x) * 180.0 / M_PI;      

    return { radius, latitude, longitude };
}
