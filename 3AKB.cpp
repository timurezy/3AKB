#include "Albert.h"
#include "Belikov.h"
#include "Cunningham.h"
#include "Simulate.h"
#include "SingleAlgorithmExecution.h"
#include "AlgorithmComparison.h"
#include "AlgorithmBenchmarking.h"
#include <array>
#include <iostream>
#include <chrono>
#include <iomanip>


int main() {
    int importedharmonics = 0;
    int nmax = 0;
    int mmax = 0;
    double radius = 6370000;
    double latitude = 0.0;
    double longitude = 0.0;
    int threads = 4;

    while (true) {

        std::string GravityModelName;
        if (selectedModel == 1) {

            GravityModelName = "EGM96";
        }
        else { GravityModelName = "EGM2008"; }

        std::cout << "\nCURRENT SETTINGS:\n"

            << "- GRAVITY MODEL " << GravityModelName << "\n"
            << "- HARMONICS: " << nmax << "\n"
            << "- COORDINATES: RADIUS = " << radius << " M, LATITUDE = " << latitude << ", LONGITUDE = " << longitude << "\n"
            << "- DOWNLOADED HARMONICS: " << importedharmonics << "\n"
            << "- THREADS = " << threads << "\n\n";

        std::cout << "SELECT METHOD:\n"
            << "1. RUN INDIVIDUAL ALGORITHM\n"
            << "2. ALGORITHM COMPARISON MODE\n"
            << "3. ALGORITHM BENCHMARKING MODE\n"
            << "4. CHANGE MAX HARMONICS\n"
            << "5. CHANGE INPUT COORDINATES\n"
            << "6. SELECT GRAVITY MODEL\n"
            << "7. IMPORT HARMONICS\n"
            << "8. NUMBER OF THREADS\n"
            << "0. EXIT\n"
            << "ENTER CHOICE: ";

        int option;
        std::cin >> option;

        if (std::cin.fail()) {
            std::cin.clear();
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            std::cout << "INVALID INPUT. TRY AGAIN.\n";
            continue;
        }
        if (option == 0) {
            freeStokes(nmax);
            std::cout << "EXIT.\n";
            break;
        }


        switch (option) {
        case 1:
            runIndividualAlgorithm(radius, latitude, longitude, nmax, mmax, threads, importedharmonics);
            break;
        case 2:
            runComparisonMode(radius, latitude, longitude, nmax, threads, importedharmonics);
            break;
        case 3:
            runBenchmarkingMode(importedharmonics);
            break;
        case 4:
            changeHarmonics(nmax, mmax, importedharmonics);
            break;
        case 5:
            changeCoordinates(radius, latitude, longitude);
            break;
        case 6:
            changeGravityModel(nmax, importedharmonics);
            break;
        case 7:
            importHarmonics(nmax, importedharmonics);
            break;
        case 8:
            changeThreads(threads);
            break;



        default:
            std::cout << "INVALID OPTION. TRY AGAIN.\n";
            break;
        }
        }

        return 0;
    }



















