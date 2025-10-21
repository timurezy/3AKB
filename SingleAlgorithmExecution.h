#pragma once
#include "Albert.h"
#include "Belikov.h"
#include "Cunningham.h"
#include "alltypes.h"
#include <array>
#include <chrono>
#include <iostream>

// Function to compute Keplerian acceleration (central body only)
std::array<double, 3> computeKeplerianAcceleration(double radius, double latitude, double longitude) {
    std::array<double, 3> Result{};
    double r = radius;
    double lat_rad = latitude * M_PI / 180.0;
    double lon_rad = longitude * M_PI / 180.0;
    
    Result[0] = -EARTH_MU * (r * cos(lat_rad) * cos(lon_rad)) / (r * r * r);
    Result[1] = -EARTH_MU * (r * cos(lat_rad) * sin(lon_rad)) / (r * r * r);
    Result[2] = -EARTH_MU * (r * sin(lat_rad)) / (r * r * r);
    
    return Result;
}

// Function to print acceleration results and timing
void printAccelerationResults(const std::array<double, 3>& Result, double total_time, const std::string& algorithm_name = "", bool is_keplerian = false) {
    std::cout << "AX = " << Result[0] << std::endl
              << "AY = " << Result[1] << std::endl
              << "AZ = " << Result[2] << std::endl
              << "TIME: " << total_time << " MS";
    
    if (is_keplerian) {
        std::cout << " (Keplerian)";
    }
    std::cout << "\n";
}

// Function to ensure harmonics are imported
void ensureHarmonicsImported(int nmax, int& importedharmonics) {
    if (importflag == 0 && nmax > 0) {
        importStokesCombined(gravityModels[selectedModel], nmax);
        std::cout << "COMBINED STOKES COEFFICIENTS IMPORTED.\n";
        importedharmonics = nmax;
    }
}

// Function to run Holmes single-thread algorithm
void runHolmesSingleThread(double radius, double latitude, double longitude, int nmax, int mmax, int& importedharmonics) {
    using namespace uniorb;
    
    auto Result = std::array<double, 3>();
    
    if (nmax == 0) {
        // Keplerian motion
        auto start = std::chrono::high_resolution_clock::now();
        Result = computeKeplerianAcceleration(radius, latitude, longitude);
        auto end = std::chrono::high_resolution_clock::now();
        double total_time = std::chrono::duration<double, std::milli>(end - start).count();
        
        printAccelerationResults(Result, total_time, "Holmes 1T", true);
        return;
    }
    
    ensureHarmonicsImported(nmax, importedharmonics);
    
    gravity_stokes GravityStokes(_c, _s, nmax, mmax, EARTH_MU, EARTH_RADIUS);
    GravityStokes.use_concurrency(1);
    
    auto start = std::chrono::high_resolution_clock::now();
    GravityStokes.get_acceleration(radius, latitude, longitude, Result);
    auto end = std::chrono::high_resolution_clock::now();
    double total_time = std::chrono::duration<double, std::milli>(end - start).count();
    
    printAccelerationResults(Result, total_time);
}

// Function to run Belikov algorithm
void runBelikov(double radius, double latitude, double longitude, int nmax, int& importedharmonics) {
    std::array<double, 3> Result{};
    
    if (nmax == 0) {
        // Keplerian motion
        auto start = std::chrono::high_resolution_clock::now();
        Result = computeKeplerianAcceleration(radius, latitude, longitude);
        auto end = std::chrono::high_resolution_clock::now();
        double total_time = std::chrono::duration<double, std::milli>(end - start).count();
        
        printAccelerationResults(Result, total_time, "Belikov", true);
        return;
    }
    
    ensureHarmonicsImported(nmax, importedharmonics);
    
    auto start = std::chrono::high_resolution_clock::now();
    gravityBelikov(radius, latitude, longitude, nmax, Result);
    auto end = std::chrono::high_resolution_clock::now();
    double total_time = std::chrono::duration<double, std::milli>(end - start).count();
    
    printAccelerationResults(Result, total_time);
}

// Function to run Cunningham algorithm
void runCunningham(double radius, double latitude, double longitude, int nmax, int& importedharmonics) {
    std::array<double, 3> Result{};
    
    if (nmax == 0) {
        // Keplerian motion
        auto start = std::chrono::high_resolution_clock::now();
        Result = computeKeplerianAcceleration(radius, latitude, longitude);
        auto end = std::chrono::high_resolution_clock::now();
        double total_time = std::chrono::duration<double, std::milli>(end - start).count();
        
        printAccelerationResults(Result, total_time, "Cunningham", true);
        return;
    }
    
    ensureHarmonicsImported(nmax, importedharmonics);
    
    auto start = std::chrono::high_resolution_clock::now();
    gravityCunningham(radius, latitude, longitude, nmax, Result);
    auto end = std::chrono::high_resolution_clock::now();
    double total_time = std::chrono::duration<double, std::milli>(end - start).count();
    
    printAccelerationResults(Result, total_time);
}

// Function to run Holmes multi-thread algorithm
void runHolmesMultiThread(double radius, double latitude, double longitude, int nmax, int mmax, int threads, int& importedharmonics) {
    using namespace uniorb;
    
    auto Result = std::array<double, 3>();
    
    if (nmax == 0) {
        // Keplerian motion
        auto start = std::chrono::high_resolution_clock::now();
        Result = computeKeplerianAcceleration(radius, latitude, longitude);
        auto end = std::chrono::high_resolution_clock::now();
        double total_time = std::chrono::duration<double, std::milli>(end - start).count();
        
        printAccelerationResults(Result, total_time, "Holmes MT", true);
        return;
    }
    
    ensureHarmonicsImported(nmax, importedharmonics);
    
    gravity_stokes GravityStokes(_c, _s, nmax, mmax, EARTH_MU, EARTH_RADIUS);
    GravityStokes.use_concurrency(threads);
    
    auto start = std::chrono::high_resolution_clock::now();
    GravityStokes.get_acceleration(radius, latitude, longitude, Result);
    auto end = std::chrono::high_resolution_clock::now();
    double total_time = std::chrono::duration<double, std::milli>(end - start).count();
    
    printAccelerationResults(Result, total_time);
}

// Function to handle harmonics change
void changeHarmonics(int& nmax, int& mmax, int& importedharmonics) {
    freeStokes(nmax);
    importedharmonics = 0;
    
    if (selectedModel == 1) {
        std::cout << "ENTER NEW NMAX (0 TO 360, 0=KEPLERIAN): ";
    }
    else {
        std::cout << "ENTER NEW NMAX (0 TO 2000, 0=KEPLERIAN): ";
    }
    
    int new_nmax;
    std::cin >> new_nmax;
    
    int max_harmonics = (selectedModel == 1) ? 360 : 2000;
    
    if (std::cin.fail() || new_nmax < 0 || new_nmax > max_harmonics) {
        std::cin.clear();
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        std::cout << "INVALID INPUT. HARMONICS NOT CHANGED.\n";
        return;
    }
    
    nmax = new_nmax;
    mmax = new_nmax;
    
    if (nmax > 0) {
        // Auto-load harmonics immediately
        importStokesCombined(gravityModels[selectedModel], nmax);
        importedharmonics = nmax;
        std::cout << "HARMONICS UPDATED AND LOADED.\n";
    }
    else {
        std::cout << "HARMONICS UPDATED (KEPLERIAN MODE).\n";
    }
}

// Function to handle coordinates change
void changeCoordinates(double& radius, double& latitude, double& longitude) {
    std::cout << "ENTER RADIUS (6000000 TO 7000000 M): ";
    double r;
    std::cin >> r;
    
    std::cout << "ENTER LATITUDE (-90 TO 90 DEG): ";
    double lat;
    std::cin >> lat;
    
    std::cout << "ENTER LONGITUDE (-180 TO 180 DEG): ";
    double lon;
    std::cin >> lon;
    
    if (std::cin.fail() || lat < -90 || lat > 90 || lon < -180 || lon > 180) {
        std::cin.clear();
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        std::cout << "INVALID INPUT. COORDINATES NOT CHANGED.\n";
        return;
    }
    
    radius = r;
    latitude = lat;
    longitude = lon;
    std::cout << "COORDINATES UPDATED.\n";
}

// Function to handle gravity model selection
void changeGravityModel(int nmax, int& importedharmonics) {
    freeStokes(nmax);
    selectGravityModel(selectedModel, gravityModels);
    importedharmonics = 0;
}

// Function to handle manual harmonics import
void importHarmonics(int nmax, int& importedharmonics) {
    freeStokes(nmax);
    importStokesCombined(gravityModels[selectedModel], nmax);
    std::cout << "COMBINED STOKES COEFFICIENTS IMPORTED.\n";
    importedharmonics = nmax;
}

// Function to handle threads change
void changeThreads(int& threads) {
    std::cout << "NUMBER OF THREADS?" << "\n";
    std::cin >> threads;
}

// Function to display individual algorithm menu
void displayIndividualAlgorithmMenu() {
    std::cout << "\nINDIVIDUAL ALGORITHM EXECUTION:\n"
              << "1. HOLMES 1T\n"
              << "2. BELIKOV\n"
              << "3. CUNNINGHAM\n"
              << "4. HOLMES MT\n"
              << "5. Back to Main Menu\n"
              << "ENTER CHOICE: ";
}

// Main individual algorithm execution handler
void runIndividualAlgorithm(double radius, double latitude, double longitude, int nmax, int mmax, 
                          int threads, int& importedharmonics) {
    bool submenu_active = true;
    
    while (submenu_active) {
        displayIndividualAlgorithmMenu();
        
        int option;
        std::cin >> option;

        if (std::cin.fail()) {
            std::cin.clear();
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            std::cout << "INVALID INPUT. TRY AGAIN.\n";
            continue;
        }

        switch (option) {
            case 1:
                runHolmesSingleThread(radius, latitude, longitude, nmax, mmax, importedharmonics);
                break;
            case 2:
                runBelikov(radius, latitude, longitude, nmax, importedharmonics);
                break;
            case 3:
                runCunningham(radius, latitude, longitude, nmax, importedharmonics);
                break;
            case 4:
                runHolmesMultiThread(radius, latitude, longitude, nmax, mmax, threads, importedharmonics);
                break;
            case 5:
                submenu_active = false;
                break;
            default:
                std::cout << "INVALID OPTION. TRY AGAIN.\n";
                break;
        }
    }
}
