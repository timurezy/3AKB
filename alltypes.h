#pragma once
#include <vector>
#include <cmath>
#include <iostream>
#define _USE_MATH_DEFINES
#include <map> 


#define EARTH_RADIUS 6378136.30 
#define EARTH_MU 398600441500000.0 

namespace Astrometric {

	struct Position {
		double x, y, z;
		Position(double _x = 0, double _y = 0, double _z = 0) : x(_x), y(_y), z(_z) {}
		std::vector<double> getCoordinates() const { return { x, y, z }; }
		double getRadius() const { return std::sqrt(x * x + y * y + z * z); }
	};

	struct StateVector {
		double x, y, z;
		StateVector(double _x = 0, double _y = 0, double _z = 0) : x(_x), y(_y), z(_z) {}
		std::vector<double> getCoordinates() const { return { x, y, z }; }
		double getRadius() const { return std::sqrt(x * x + y * y + z * z); }
	};

	struct StateVectorSph {
		double r, t, u; 

		StateVectorSph(double _r = 0, double _t = 0, double _u = 0)
			: r(_r), t(_t), u(_u) {
		}

		std::vector<double> getCoordinates() const {
			return { r, t, u };
		}

		double getRadius() const {
			return r;
		}
	};


	struct Matrix {
		Position pos;
		Matrix(const Position& _p) : pos(_p) {}
	};

	namespace Text {
		inline void error(const char* msg) {
			std::cerr << msg << std::endl;
		}
	}

	void gravityBelikov(double r, double lat, double lon, int n, std::array<double, 3>& Result);
	void gravityCunningham(double r, double lat, double lon, int n, std::array<double, 3>& Result);

}

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
