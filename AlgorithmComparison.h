#pragma once
#include "Albert.h"
#include "Belikov.h"
#include "Cunningham.h"
#include "Simulate.h"
#include "alltypes.h"
#include <array>
#include <chrono>
#include <iostream>
#include <fstream>
#include <iomanip>

// Function to display the comparison mode submenu
void displayComparisonMenu(const std::string& gravityModelName, int nmax, double radius, 
                          double latitude, double longitude, int importedharmonics, 
                          int num_runs, int threads) {
    std::cout << "\nCURRENT SETTINGS:\n"
              << "- GRAVITY MODEL " << gravityModelName << "\n"
              << "- HARMONICS: " << nmax << "\n"
              << "- COORDINATES: RADIUS = " << radius << " M, LATITUDE = " << latitude 
              << ", LONGITUDE = " << longitude << "\n"
              << "- DOWNLOADED HARMONICS: " << importedharmonics << "\n"
              << "- NUMBER OF RUNS: " << num_runs << "\n"
              << "- THREADS = " << threads << "\n\n";

    std::cout << "\nCOMPARISON MODE MENU:\n"
              << "1. Random Data Comparison\n"
              << "2. Sequential Data Comparison\n"
              << "3. Change Number of Runs\n"
              << "4. Back to Main Menu\n"
              << "Enter choice: ";
}



// Function to handle random data comparison with file output
void runRandomDataComparison(int nmax, int num_runs, int threads, int& importedharmonics) {
    std::cout << "Random Data Mode to file (formatted output) selected.\n";

    if (importflag == 0) {
        importStokesCombined(gravityModels[selectedModel], nmax);
        importedharmonics = nmax;
    }

    std::ofstream file("results_RANDOM.csv");
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file for writing.\n";
        return;
    }

    file << "Algorithm,Run,Time (ms),R,Latitude,Longitude,ax,ay,az\n";

    double total_time_belikov = 0.0;
    double total_time_cunningham = 0.0;
    double total_time_stokes_single = 0.0;
    double total_time_stokes_multi = 0.0;

    std::array<double, 3> resB = {0.0, 0.0, 0.0};
    std::array<double, 3> resC = {0.0, 0.0, 0.0};
    std::array<double, 3> resS1 = {0.0, 0.0, 0.0};
    std::array<double, 3> resSM = {0.0, 0.0, 0.0};

    std::cout << "Algorithm is running...\n";

    for (int run = 1; run <= num_runs; ++run) {
        auto coords = generateRandomCoordinates();
        double radius = coords[0];
        double latitude = coords[1];
        double longitude = coords[2];

        // === Belikov ===
        auto startB = std::chrono::high_resolution_clock::now();
        gravityBelikov(radius, latitude, longitude, nmax, resB);
        simulate_integrator_and_sofa(27);
        auto endB = std::chrono::high_resolution_clock::now();
        double timeB = std::chrono::duration<double, std::milli>(endB - startB).count();
        file << "Belikov," << run << "," << timeB << ","
             << radius << "," << latitude << "," << longitude << ","
             << resB[0] << "," << resB[1] << "," << resB[2] << "\n";
        total_time_belikov += timeB;

        // === Cunningham ===
        auto startC = std::chrono::high_resolution_clock::now();
        gravityCunningham(radius, latitude, longitude, nmax, resC);
        simulate_integrator_and_sofa(27);
        auto endC = std::chrono::high_resolution_clock::now();
        double timeC = std::chrono::duration<double, std::milli>(endC - startC).count();
        file << "Cunningham," << run << "," << timeC << ","
             << radius << "," << latitude << "," << longitude << ","
             << resC[0] << "," << resC[1] << "," << resC[2] << "\n";
        total_time_cunningham += timeC;

        // === Stokes (1-thread) ===
        using namespace uniorb;
        gravity_stokes GravityStokes(_c, _s, nmax, nmax, EARTH_MU, EARTH_RADIUS);
        GravityStokes.use_concurrency(1);

        auto startS1 = std::chrono::high_resolution_clock::now();
        GravityStokes.get_acceleration(radius, latitude, longitude, resS1);
        simulate_integrator_and_sofa(27);
        auto endS1 = std::chrono::high_resolution_clock::now();
        double timeS1 = std::chrono::duration<double, std::milli>(endS1 - startS1).count();
        file << "Stokes_1thread," << run << "," << timeS1 << ","
             << radius << "," << latitude << "," << longitude << ","
             << resS1[0] << "," << resS1[1] << "," << resS1[2] << "\n";
        total_time_stokes_single += timeS1;

        // === Stokes (multi-thread) ===
        GravityStokes.use_concurrency(threads);

        auto startSM = std::chrono::high_resolution_clock::now();
        GravityStokes.get_acceleration(radius, latitude, longitude, resSM);
        simulate_integrator_and_sofa(27);
        auto endSM = std::chrono::high_resolution_clock::now();
        double timeSM = std::chrono::duration<double, std::milli>(endSM - startSM).count();
        file << "Stokes_multithread," << run << "," << timeSM << ","
             << radius << "," << latitude << "," << longitude << ","
             << resSM[0] << "," << resSM[1] << "," << resSM[2] << "\n";
        total_time_stokes_multi += timeSM;

        int progress = static_cast<int>(100.0 * run / num_runs);
        if (progress % 5 == 0)
            std::cout << "\rProgress: " << progress << "%" << std::flush;
    }

    std::cout << "\rProgress: 100%\nComputation complete.\n";

    double speedup = total_time_stokes_single / total_time_stokes_multi;
    double efficiency = speedup / threads;

    file << "\nTotal time (ms):\n";
    file << "Belikov," << total_time_belikov << "\n";
    file << "Cunningham," << total_time_cunningham << "\n";
    file << "Stokes_1thread," << total_time_stokes_single << "\n";
    file << "Stokes_multithread," << total_time_stokes_multi << "\n";

    file << "Speedup," << speedup << "\n";
    file << "Efficiency," << efficiency << "\n";

    //file << "\nFinal coordinates and accelerations:\n";
    //file << "R=" << coords[0] << " LAT=" << coords[1] << " LON=" << coords[2] << "\n";

    file << std::fixed << std::setprecision(10);
    file << "Final accelerations (m/s^2):\n";
    file << "Belikov,"
         << resB[0] << "," << resB[1] << "," << resB[2] << "\n";
    file << "Cunningham,"
         << resC[0] << "," << resC[1] << "," << resC[2] << "\n";
    file << "Holmes_1thread,"
         << resS1[0] << "," << resS1[1] << "," << resS1[2] << "\n";
    file << "Holmes_multithread,"
         << resSM[0] << "," << resSM[1] << "," << resSM[2] << "\n";

    file.close();
    std::cout << "Results saved to results_RANDOM.csv\n";
}

// Function to handle sequential data comparison with file output
void runSequentialDataComparison(double radius, double latitude, double longitude, int nmax, 
                                int num_runs, int threads, int& importedharmonics) {
    std::cout << "Sequential Data Mode to file (formatted output) selected.\n";

    if (importflag == 0) {
        importStokesCombined(gravityModels[selectedModel], nmax);
        importedharmonics = nmax;
    }

    std::ofstream file("results_SEQUENTIAL.csv");
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file for writing.\n";
        return;
    }

    file << "Algorithm,Run,Time (ms),R,Latitude,Longitude,ax,ay,az\n";

    std::array<double, 3> xyz = RLatLonToXYZ(radius, latitude, longitude);
    std::array<double, 3> v = random_velocity(radius);

    std::array<double, 3> sphericalB = {0,0,0};
    std::array<double, 3> sphericalC = {0,0,0};
    std::array<double, 3> sphericalH1 = {0,0,0};
    std::array<double, 3> sphericalHm = {0,0,0};

    double fulltimeB = 0, fulltimeC = 0, fulltimeH1 = 0, fulltimeHm = 0;
    double s = 0, e = 0;

    double xB = xyz[0], xC = xyz[0], xH1 = xyz[0], xHm = xyz[0];
    double yB = xyz[1], yC = xyz[1], yH1 = xyz[1], yHm = xyz[1];
    double zB = xyz[2], zC = xyz[2], zH1 = xyz[2], zHm = xyz[2];

    double vxB = v[0], vxC = v[0], vxH1 = v[0], vxHm = v[0];
    double vyB = v[1], vyC = v[1], vyH1 = v[1], vyHm = v[1];
    double vzB = v[2], vzC = v[2], vzH1 = v[2], vzHm = v[2];

    double radiusB = radius, radiusC = radius, radiusH1 = radius, radiusHm = radius;
    double latitudeB = latitude, latitudeC = latitude, latitudeH1 = latitude, latitudeHm = latitude;
    double longitudeB = longitude, longitudeC = longitude, longitudeH1 = longitude, longitudeHm = longitude;

    double deltat = 0.001;

    std::cout << "Algorithm is running...\n";

    for (int run = 1; run <= num_runs; ++run) {
        // === Belikov ===
        std::array<double, 3> ResultBelikov{};
        auto startB = std::chrono::high_resolution_clock::now();
        gravityBelikov(radiusB, latitudeB, longitudeB, nmax, ResultBelikov);
        simulate_integrator_and_sofa(27);
        auto endB = std::chrono::high_resolution_clock::now();
        double timeB = std::chrono::duration<double, std::milli>(endB - startB).count();
        fulltimeB += timeB;

        file << "Belikov," << run << "," << timeB << ","
             << radiusB << "," << latitudeB << "," << longitudeB << ","
             << ResultBelikov[0] << "," << ResultBelikov[1] << "," << ResultBelikov[2] << "\n";

        // Integration
        xB += vxB * deltat + 0.5 * ResultBelikov[0] * deltat * deltat;
        yB += vyB * deltat + 0.5 * ResultBelikov[1] * deltat * deltat;
        zB += vzB * deltat + 0.5 * ResultBelikov[2] * deltat * deltat;
        vxB += ResultBelikov[0] * deltat;
        vyB += ResultBelikov[1] * deltat;
        vzB += ResultBelikov[2] * deltat;
        sphericalB = XYZtoRLatLon(xB, yB, zB);
        radiusB = sphericalB[0];
        latitudeB = sphericalB[1];
        longitudeB = sphericalB[2];

        // === Cunningham ===
        std::array<double, 3> ResultCunningham{};
        auto startC = std::chrono::high_resolution_clock::now();
        gravityCunningham(radiusC, latitudeC, longitudeC, nmax, ResultCunningham);
        simulate_integrator_and_sofa(27);
        auto endC = std::chrono::high_resolution_clock::now();
        double timeC = std::chrono::duration<double, std::milli>(endC - startC).count();
        fulltimeC += timeC;

        file << "Cunningham," << run << "," << timeC << ","
             << radiusC << "," << latitudeC << "," << longitudeC << ","
             << ResultCunningham[0] << "," << ResultCunningham[1] << "," << ResultCunningham[2] << "\n";

        xC += vxC * deltat + 0.5 * ResultCunningham[0] * deltat * deltat;
        yC += vyC * deltat + 0.5 * ResultCunningham[1] * deltat * deltat;
        zC += vzC * deltat + 0.5 * ResultCunningham[2] * deltat * deltat;
        vxC += ResultCunningham[0] * deltat;
        vyC += ResultCunningham[1] * deltat;
        vzC += ResultCunningham[2] * deltat;
        sphericalC = XYZtoRLatLon(xC, yC, zC);
        radiusC = sphericalC[0];
        latitudeC = sphericalC[1];
        longitudeC = sphericalC[2];

        // === Stokes (1-thread) ===
        using namespace uniorb;
        std::array<double, 3> ResultStokesONE{};
        gravity_stokes GravityStokes(_c, _s, nmax, nmax, EARTH_MU, EARTH_RADIUS);
        GravityStokes.use_concurrency(1);

        auto startS = std::chrono::high_resolution_clock::now();
        GravityStokes.get_acceleration(radiusH1, latitudeH1, longitudeH1, ResultStokesONE);
        simulate_integrator_and_sofa(27);
        auto endS = std::chrono::high_resolution_clock::now();
        double timeS = std::chrono::duration<double, std::milli>(endS - startS).count();
        fulltimeH1 += timeS;

        file << "Stokes_1thread," << run << "," << timeS << ","
             << radiusH1 << "," << latitudeH1 << "," << longitudeH1 << ","
             << ResultStokesONE[0] << "," << ResultStokesONE[1] << "," << ResultStokesONE[2] << "\n";

        xH1 += vxH1 * deltat + 0.5 * ResultStokesONE[0] * deltat * deltat;
        yH1 += vyH1 * deltat + 0.5 * ResultStokesONE[1] * deltat * deltat;
        zH1 += vzH1 * deltat + 0.5 * ResultStokesONE[2] * deltat * deltat;
        vxH1 += ResultStokesONE[0] * deltat;
        vyH1 += ResultStokesONE[1] * deltat;
        vzH1 += ResultStokesONE[2] * deltat;
        sphericalH1 = XYZtoRLatLon(xH1, yH1, zH1);
        radiusH1 = sphericalH1[0];
        latitudeH1 = sphericalH1[1];
        longitudeH1 = sphericalH1[2];

        // === Stokes (multi-thread) ===
        std::array<double, 3> ResultStokesMULTI{};
        GravityStokes.use_concurrency(threads);

        auto startSM = std::chrono::high_resolution_clock::now();
        GravityStokes.get_acceleration(radiusHm, latitudeHm, longitudeHm, ResultStokesMULTI);
        simulate_integrator_and_sofa(27);
        auto endSM = std::chrono::high_resolution_clock::now();
        double timeSM = std::chrono::duration<double, std::milli>(endSM - startSM).count();
        fulltimeHm += timeSM;

        file << "Stokes_multithread," << run << "," << timeSM << ","
             << radiusHm << "," << latitudeHm << "," << longitudeHm << ","
             << ResultStokesMULTI[0] << "," << ResultStokesMULTI[1] << "," << ResultStokesMULTI[2] << "\n";

        xHm += vxHm * deltat + 0.5 * ResultStokesMULTI[0] * deltat * deltat;
        yHm += vyHm * deltat + 0.5 * ResultStokesMULTI[1] * deltat * deltat;
        zHm += vzHm * deltat + 0.5 * ResultStokesMULTI[2] * deltat * deltat;
        vxHm += ResultStokesMULTI[0] * deltat;
        vyHm += ResultStokesMULTI[1] * deltat;
        vzHm += ResultStokesMULTI[2] * deltat;
        sphericalHm = XYZtoRLatLon(xHm, yHm, zHm);
        radiusHm = sphericalHm[0];
        latitudeHm = sphericalHm[1];
        longitudeHm = sphericalHm[2];

        int progress = static_cast<int>(100.0 * run / num_runs);
        if (progress % 5 == 0)
            std::cout << "\rProgress: " << progress << "%" << std::flush;
    }

    std::cout << "\rProgress: 100%\nComputation complete.\n";

    s = fulltimeH1 / fulltimeHm;
    e = s / threads;

    file << "\nTotal time (ms):\n";
    file << "Belikov," << fulltimeB << "\n";
    file << "Cunningham," << fulltimeC << "\n";
    file << "Stokes_1thread," << fulltimeH1 << "\n";
    file << "Stokes_multithread," << fulltimeHm << "\n";
    file << "Speedup," << s << "\n";
    file << "Efficiency," << e << "\n";

    file << "\nFinal coordinates (R, LAT, LON):\n";
    file << "Belikov," << radiusB << "," << latitudeB << "," << longitudeB << "\n";
    file << "Cunningham," << radiusC << "," << latitudeC << "," << longitudeC << "\n";
    file << "Stokes_1thread," << radiusH1 << "," << latitudeH1 << "," << longitudeH1 << "\n";
    file << "Stokes_multithread," << radiusHm << "," << latitudeHm << "," << longitudeHm << "\n";

    file.close();
    std::cout << "Results saved to results_SEQUENTIAL.csv\n";
}

// Function to handle changing number of runs
bool changeNumberOfRuns(int& num_runs) {
    std::cout << "Enter new number of runs: ";
    std::cin >> num_runs;
    if (std::cin.fail() || num_runs < 1) {
        std::cin.clear();
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        std::cout << "INVALID INPUT. NUMBER OF RUNS NOT CHANGED.\n";
        return false;
    }
    else {
        std::cout << "Number of runs updated to " << num_runs << ".\n";
        return true;
    }
}

// Main comparison mode handler function
void runComparisonMode(double& radius, double& latitude, double& longitude, int nmax, 
                      int threads, int& importedharmonics) {
    std::string gravityModelName = (selectedModel == 1) ? "EGM96" : "EGM2008";
    bool submenu_active = true;
    int num_runs = 100;
    
    while (submenu_active) {
        displayComparisonMenu(gravityModelName, nmax, radius, latitude, longitude, 
                            importedharmonics, num_runs, threads);
        
        int sub_option;
        std::cin >> sub_option;

        if (std::cin.fail()) {
            std::cin.clear();
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            std::cout << "INVALID INPUT. TRY AGAIN.\n";
            continue;
        }

        switch (sub_option) {
            case 1:
                runRandomDataComparison(nmax, num_runs, threads, importedharmonics);
                break;
            case 2:
                runSequentialDataComparison(radius, latitude, longitude, nmax, num_runs, threads, importedharmonics);
                break;
            case 3:
                changeNumberOfRuns(num_runs);
                break;
            case 4:
                submenu_active = false;
                break;
            default:
                std::cout << "INVALID OPTION. TRY AGAIN.\n";
                break;
        }
    }
}