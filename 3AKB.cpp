#include "Albert.h"
#include "Belikov.h"
#include "Cunningham.h"
#include <array>
#include <iostream>
#include <chrono>

int main() {
    int nmax = 16;
    int mmax = 16;
    double radius = 6370000;
    double latitude = 0.0;
    double longitude = 0.0;

    while (true) {

        std::string GravityModelName;
        if (selectedModel == 1) {

            GravityModelName = "EGM96";
        }
        else { GravityModelName = "EGM2008"; }

        std::cout << "\nCURRENT SETTINGS:\n"

            << "- GRAVITY MODEL " << GravityModelName << "\n"
            << "- HARMONICS: " << nmax << "\n"
            << "- COORDINATES: RADIUS = " << radius << " M, LATITUDE = " << latitude << ", LONGITUDE = " << longitude << "\n\n";

        std::cout << "SELECT METHOD:\n"
            << "1. HOLMES 1T\n"
            << "2. BELIKOV\n"
            << "3. CUNNINGHAM\n"
            << "4. HOLMES MT\n"
            << "5. CHANGE MAX HARMONICS\n"
            << "6. CHANGE INPUT COORDINATES\n"
            << "7. SELECT GRAVITY MODEL\n"
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
            std::cout << "EXIT.\n";
            break;
        }

        switch (option) {
        case 1: {
            using namespace uniorb;

            auto Result = std::array<double, 3>();
            auto GravityStokes = gravity_stokes();
            GravityStokes.import(gravityModels[selectedModel], nmax, mmax);
            GravityStokes.use_concurrency(1);

            auto start = std::chrono::high_resolution_clock::now();
            GravityStokes.get_acceleration(radius, latitude, longitude, Result);
            auto end = std::chrono::high_resolution_clock::now();
            double total_time = std::chrono::duration<double, std::milli>(end - start).count();  // ✅ изменено

            std::cout
                << "AX = " << Result[0] << std::endl
                << "AY = " << Result[1] << std::endl
                << "AZ = " << Result[2] << std::endl
                << "TIME: " << total_time << " MS\n";  // ✅ изменено
            break;
        }
        case 2: {

            std::array<double, 3> Result{};
            importStokesBelikov(gravityModels[selectedModel], nmax);

            auto start = std::chrono::high_resolution_clock::now();
            gravityBelikov(radius, latitude, longitude, nmax, Result);
            auto end = std::chrono::high_resolution_clock::now();
            double total_time = std::chrono::duration<double, std::milli>(end - start).count();  // ✅ изменено

            std::cout
                << "AX = " << Result[0] << "\n"
                << "AY = " << Result[1] << "\n"
                << "AZ = " << Result[2] << "\n"
                << "TIME: " << total_time << " MS\n";  // ✅ изменено
            break;
        }
        case 3: {

            std::array<double, 3> Result{};
            importStokesCunningham(gravityModels[selectedModel], nmax);

            auto start = std::chrono::high_resolution_clock::now();
            gravityCunningham(radius, latitude, longitude, nmax, Result);
            auto end = std::chrono::high_resolution_clock::now();
            double total_time = std::chrono::duration<double, std::milli>(end - start).count();  // ✅ изменено

            std::cout
                << "AX = " << Result[0] << "\n"
                << "AY = " << Result[1] << "\n"
                << "AZ = " << Result[2] << "\n"
                << "TIME: " << total_time << " MS\n";  // ✅ изменено
            break;
        }
        case 4: {
            using namespace uniorb;

            auto Result = std::array<double, 3>();
            auto GravityStokes = gravity_stokes();
            GravityStokes.import(gravityModels[selectedModel], nmax, mmax);
            GravityStokes.use_concurrency(4);

            auto start = std::chrono::high_resolution_clock::now();
            GravityStokes.get_acceleration(radius, latitude, longitude, Result);
            auto end = std::chrono::high_resolution_clock::now();
            double total_time = std::chrono::duration<double, std::milli>(end - start).count();  // ✅ изменено

            std::cout
                << "AX = " << Result[0] << std::endl
                << "AY = " << Result[1] << std::endl
                << "AZ = " << Result[2] << std::endl
                << "TIME: " << total_time << " MS\n";  // ✅ изменено
            break;
        }
        case 5: {
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
            selectGravityModel(selectedModel, gravityModels);
            break;
        }
        default:
            std::cout << "INVALID OPTION. TRY AGAIN.\n";
            break;
        }
    }

    return 0;
}


















