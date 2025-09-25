#include "Albert.h"
#include "Belikov.h"
#include "Cunningham.h"
#include "Simulate.h"
#include <array>
#include <iostream>
#include <chrono>



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
            << "1. HOLMES 1T\n"
            << "2. BELIKOV\n"
            << "3. CUNNINGHAM\n"
            << "4. HOLMES MT\n"
            << "5. CHANGE MAX HARMONICS\n"
            << "6. CHANGE INPUT COORDINATES\n"
            << "7. SELECT GRAVITY MODEL\n"
            << "8. IMPORT HARMONICS\n"
            << "9. NUMBER OF THREADS\n"
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
        case 1: {
            using namespace uniorb;

            if (importflag == 0) {
                importStokesCombined(gravityModels[selectedModel], nmax);
                std::cout << "COMBINED STOKES COEFFICIENTS IMPORTED.\n";
                importedharmonics = nmax;
            }
            auto Result = std::array<double, 3>();
            //auto GravityStokes = gravity_stokes();
            gravity_stokes GravityStokes(_c, _s, nmax, mmax, EARTH_MU, EARTH_RADIUS);

            GravityStokes.use_concurrency(1);

            auto start = std::chrono::high_resolution_clock::now();
            GravityStokes.get_acceleration(radius, latitude, longitude, Result);
            auto end = std::chrono::high_resolution_clock::now();
            double total_time = std::chrono::duration<double, std::milli>(end - start).count();

            std::cout
                << "AX = " << Result[0] << std::endl
                << "AY = " << Result[1] << std::endl
                << "AZ = " << Result[2] << std::endl
                << "TIME: " << total_time << " MS\n";
            break;
        }
        case 2: {

            std::array<double, 3> Result{};
            if (importflag == 0) {
                importStokesCombined(gravityModels[selectedModel], nmax);
                std::cout << "COMBINED STOKES COEFFICIENTS IMPORTED.\n";
                importedharmonics = nmax;
            }
            auto start = std::chrono::high_resolution_clock::now();
            gravityBelikov(radius, latitude, longitude, nmax, Result);
            auto end = std::chrono::high_resolution_clock::now();
            double total_time = std::chrono::duration<double, std::milli>(end - start).count();

            std::cout
                << "AX = " << Result[0] << "\n"
                << "AY = " << Result[1] << "\n"
                << "AZ = " << Result[2] << "\n"
                << "TIME: " << total_time << " MS\n";
            break;
        }
        case 3: {

            std::array<double, 3> Result{};
            if (importflag == 0) {
                importStokesCombined(gravityModels[selectedModel], nmax);
                std::cout << "COMBINED STOKES COEFFICIENTS IMPORTED.\n";
                importedharmonics = nmax;
            }
            auto start = std::chrono::high_resolution_clock::now();
            gravityCunningham(radius, latitude, longitude, nmax, Result);
            auto end = std::chrono::high_resolution_clock::now();
            double total_time = std::chrono::duration<double, std::milli>(end - start).count();

            std::cout
                << "AX = " << Result[0] << "\n"
                << "AY = " << Result[1] << "\n"
                << "AZ = " << Result[2] << "\n"
                << "TIME: " << total_time << " MS\n";
            break;
        }
        case 4: {

            using namespace uniorb;

            auto Result = std::array<double, 3>();
            gravity_stokes GravityStokes(_c, _s, nmax, mmax, EARTH_MU, EARTH_RADIUS);
            if (importflag == 0) {
                importStokesCombined(gravityModels[selectedModel], nmax);
                std::cout << "COMBINED STOKES COEFFICIENTS IMPORTED.\n";
                importedharmonics = nmax;
            }
            GravityStokes.use_concurrency(threads);

            auto start = std::chrono::high_resolution_clock::now();
            GravityStokes.get_acceleration(radius, latitude, longitude, Result);
            auto end = std::chrono::high_resolution_clock::now();
            double total_time = std::chrono::duration<double, std::milli>(end - start).count();

            std::cout
                << "AX = " << Result[0] << std::endl
                << "AY = " << Result[1] << std::endl
                << "AZ = " << Result[2] << std::endl
                << "TIME: " << total_time << " MS\n";
            break;
        }
        case 5: {
            freeStokes(nmax);
            importedharmonics = 0;
            if (selectedModel == 1) {
                std::cout << "ENTER NEW NMAX (1 TO 360): ";
            }
            else {
                std::cout << "ENTER NEW NMAX (1 TO 2000): ";
            }
            int new_nmax;
            std::cin >> new_nmax;
            if (selectedModel == 1) {
                if (std::cin.fail() || new_nmax < 1 || new_nmax > 360) {
                    std::cin.clear();
                    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                    std::cout << "INVALID INPUT. HARMONICS NOT CHANGED.\n";
                }
                else {
                    nmax = new_nmax;
                    mmax = new_nmax;
                    std::cout << "HARMONICS UPDATED.\n";
                }
            }
            else {
                if (std::cin.fail() || new_nmax < 1 || new_nmax > 2000) {
                    std::cin.clear();
                    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                    std::cout << "INVALID INPUT. HARMONICS NOT CHANGED.\n";
                }
                else {
                    nmax = new_nmax;
                    mmax = new_nmax;
                    std::cout << "HARMONICS UPDATED.\n";
                }
            }
            break;
        }
        case 6: {
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
            }
            else {
                radius = r;
                latitude = lat;
                longitude = lon;
                std::cout << "COORDINATES UPDATED.\n";
            }
            break;
        }
        case 7: {
            freeStokes(nmax);
            selectGravityModel(selectedModel, gravityModels);
            importedharmonics = 0;
            break;
        }
        case 8: {
            freeStokes(nmax);
            importStokesCombined(gravityModels[selectedModel], nmax);
            std::cout << "COMBINED STOKES COEFFICIENTS IMPORTED.\n";
            importedharmonics = nmax;
            break;
        }
        case 9: {
            std::cout << "NUMBER OF THREADS?" << "\n";
            std::cin >> threads;
            break;
        }

        case 10: {
            bool submenu_active = true;
            int num_runs = 100;
            while (submenu_active) {




                std::cout << "\nCURRENT SETTINGS:\n"

                    << "- GRAVITY MODEL " << GravityModelName << "\n"
                    << "- HARMONICS: " << nmax << "\n"
                    << "- COORDINATES: RADIUS = " << radius << " M, LATITUDE = " << latitude << ", LONGITUDE = " << longitude << "\n"
                    << "- DOWNLOADED HARMONICS: " << importedharmonics << "\n"
                    << "- THREADS = " << threads << "\n\n";




                std::cout << "\nDATA MODE MENU:\n"
                    << "1. Random Data Mode\n"
                    << "2. Sequential Data Mode\n"
                    << "3. Change Number of Runs\n"
                    << "4. Back to Main Menu\n"
                    << "Enter choice: ";
                int sub_option;
                std::cin >> sub_option;

                if (std::cin.fail()) {
                    std::cin.clear();
                    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                    std::cout << "INVALID INPUT. TRY AGAIN.\n";
                    continue;
                }

                switch (sub_option) {
                case 1: {
                    std::cout << "Random Data Mode selected.\n";


                    if (importflag == 0) {
                        importStokesCombined(gravityModels[selectedModel], nmax);
                        std::cout << "COMBINED STOKES COEFFICIENTS IMPORTED.\n";
                        importedharmonics = nmax;
                    }


                    for (int i = 0; i < num_runs; i++) {
                        auto coords = generateRandomCoordinates();
                        radius = coords[0];
                        latitude = coords[1];
                        longitude = coords[2];

                        std::cout << "\nRun " << i + 1
                            << ": R=" << radius
                            << " LAT=" << latitude
                            << " LON=" << longitude << "\n";

                        // === Belikov ===
                        std::array<double, 3> ResultBelikov{};
                        auto startB = std::chrono::high_resolution_clock::now();
                        gravityBelikov(radius, latitude, longitude, nmax, ResultBelikov);
                        auto endB = std::chrono::high_resolution_clock::now();
                        double timeB = std::chrono::duration<double, std::milli>(endB - startB).count();

                        std::cout << "Belikov: AX=" << ResultBelikov[0]
                            << " AY=" << ResultBelikov[1]
                            << " AZ=" << ResultBelikov[2]
                            << " TIME=" << timeB << " ms\n";

                        // === Cunningham ===
                        std::array<double, 3> ResultCunningham{};
                        auto startC = std::chrono::high_resolution_clock::now();
                        gravityCunningham(radius, latitude, longitude, nmax, ResultCunningham);
                        auto endC = std::chrono::high_resolution_clock::now();
                        double timeC = std::chrono::duration<double, std::milli>(endC - startC).count();

                        std::cout << "Cunningham: AX=" << ResultCunningham[0]
                            << " AY=" << ResultCunningham[1]
                            << " AZ=" << ResultCunningham[2]
                            << " TIME=" << timeC << " ms\n";

                        // === Stokes (Holmes 1-thread) ===
                        using namespace uniorb;
                        std::array<double, 3> ResultStokesONE{};
                        gravity_stokes GravityStokes(_c, _s, nmax, mmax, EARTH_MU, EARTH_RADIUS);
                        GravityStokes.use_concurrency(1);

                        auto startS = std::chrono::high_resolution_clock::now();
                        GravityStokes.get_acceleration(radius, latitude, longitude, ResultStokesONE);
                        auto endS = std::chrono::high_resolution_clock::now();
                        double timeS = std::chrono::duration<double, std::milli>(endS - startS).count();

                        std::cout << "Stokes: AX=" << ResultStokesONE[0]
                            << " AY=" << ResultStokesONE[1]
                            << " AZ=" << ResultStokesONE[2]
                            << " TIME=" << timeS << " ms\n";

                        std::array<double, 3> ResultStokesMULTI{};

                        GravityStokes.use_concurrency(threads);

                        auto startSM = std::chrono::high_resolution_clock::now();
                        GravityStokes.get_acceleration(radius, latitude, longitude, ResultStokesMULTI);
                        auto endSM = std::chrono::high_resolution_clock::now();
                        double timeSM = std::chrono::duration<double, std::milli>(endSM - startSM).count();

                        std::cout << "Stokes: AX=" << ResultStokesMULTI[0]
                            << " AY=" << ResultStokesMULTI[1]
                            << " AZ=" << ResultStokesMULTI[2]
                            << " TIME=" << timeSM << " ms\n";

                    }
                    break;
                }

                case 2:
                    std::cout << "Sequential Data Mode selected.\n";
                    
                    break;
                case 3:
                    std::cout << "Enter new number of runs: ";
                    std::cin >> num_runs;
                    if (std::cin.fail() || num_runs < 1) {
                        std::cin.clear();
                        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                        std::cout << "INVALID INPUT. NUMBER OF RUNS NOT CHANGED.\n";
                    }
                    else {
                        std::cout << "Number of runs updated to " << num_runs << ".\n";
                    }
                    break;
                case 4:
                    submenu_active = false; 
                    break;
                default:
                    std::cout << "INVALID OPTION. TRY AGAIN.\n";
                    break;
                }
                }
                break;
            }



        default:
            std::cout << "INVALID OPTION. TRY AGAIN.\n";
            break;
        }
        }

        return 0;
    }



















