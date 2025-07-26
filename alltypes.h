#pragma once
#include <vector>
#include <cmath>
#include <iostream>
#define _USE_MATH_DEFINES
#include <map> 


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
