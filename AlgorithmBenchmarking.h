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
#include <thread>
#include <vector>
#include <string>


std::string cleanName(std::string name) {
    // Çàìåíÿåì ïðîáåëû íà '_'
    std::replace(name.begin(), name.end(), ' ', '_');
    // Óäàëÿåì ïåðåíîñû ñòðîê è âîçâðàòû êàðåòêè
    name.erase(std::remove_if(name.begin(), name.end(), [](char c) { return c == '\n' || c == '\r'; }), name.end());
    // Óäàëÿåì ìíîæåñòâåííûå '_': çàìåíÿåì '__' íà '_', ïîêà îíè åñòü
    size_t pos;
    while ((pos = name.find("__")) != std::string::npos) {
        name.replace(pos, 2, "_");
    }
    // Óäàëÿåì âåäóùèå è trailing '_'
    if (!name.empty() && name[0] == '_') name.erase(0, 1);
    if (!name.empty() && name.back() == '_') name.pop_back();
    return name;
}

// Ôóíêöèÿ äëÿ ïîëó÷åíèÿ èìåíè ìîäåëè ïðîöåññîðà
std::string getProcessorName() {
    std::string processorName;
#ifdef _WIN32
    // Âûçîâ wmic äëÿ Windows
    std::stringstream cmd;
    cmd << "wmic cpu get name /value";
    FILE* pipe = _popen(cmd.str().c_str(), "r");
    if (pipe) {
        char buffer[256];
        std::string result;
        while (fgets(buffer, sizeof(buffer), pipe) != NULL) {
            result += buffer;
        }
        _pclose(pipe);
        // Ïàðñèíã: èùåì "Name=" è èçâëåêàåì çíà÷åíèå
        size_t pos = result.find("Name=");
        if (pos != std::string::npos) {
            pos += 5; // Ïðîïóñòèòü "Name="
            size_t end = result.find("\n", pos);
            processorName = result.substr(pos, end - pos);
        }
    }
#else
    // Äëÿ Linux: èñïîëüçóéòå popen
    FILE* pipe = popen("grep 'model name' /proc/cpuinfo | head -1 | cut -d: -f2 | tr -d '\\n'", "r");
    if (pipe) {
        char buffer[256];
        if (fgets(buffer, sizeof(buffer), pipe) != NULL) {
            processorName = std::string(buffer);
        }
        pclose(pipe);
    }
#endif
    if (processorName.empty()) {
        processorName = "Unknown_CPU";
    }
    // Î÷èùàåì èìÿ
    return cleanName(processorName);
}


// Get the number of hardware threads available on the system
int getMaxThreads() {
    return std::thread::hardware_concurrency();
}

// Function to run benchmarks across all harmonics and thread counts
void runComprehensiveBenchmark(int num_runs, int& importedharmonics) {

    bool isUseImitation = false;

    std::cout << "=== COMPREHENSIVE ALGORITHM BENCHMARK ===\n";
    
    int maxThreads = getMaxThreads();
    std::cout << "Detected " << maxThreads << " hardware threads.\n";
    
    int maxHarmonics = (selectedModel == 1) ? 360 : 2000;
    std::string gravityModelName = (selectedModel == 1) ? "EGM96" : "EGM2008";
    
    std::cout << "Selected model: " << gravityModelName << "\n";
    std::cout << "Max harmonics: " << maxHarmonics << "\n";
    std::cout << "Number of runs per test: " << num_runs << "\n";
    std::cout << "This will take a considerable amount of time...\n\n";
    
    std::cout << "Do you want to proceed? (y/n): ";
    char confirm;
    std::cin >> confirm;
    if (confirm != 'y' && confirm != 'Y') {
        std::cout << "Benchmark cancelled.\n";
        return;
    }
    
    // Create output file

    // Ïîëó÷àåì èìÿ ïðîöåññîðà
    std::string processorName = getProcessorName();

    // Ãåíåðèðóåì èìÿ ôàéëà
    std::stringstream filename;
    filename << "benchmark_" << processorName << "_" << maxThreads << "_" << num_runs << ".csv";

    // Îòêðûâàåì ôàéë
    std::ofstream file(filename.str());

    //std::ofstream file("benchmark_comprehensive.csv");
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file for writing.\n";
        return;
    }
    
    // Write CSV header
    file << "Algorithm,Threads,Harmonics,AvgTime_ms,StdDev_ms,TotalTime_ms,Speedup,Efficiency\n";
    
    std::cout << "\nStarting benchmark...\n";
    
    // Iterate through harmonics from 2 to maxHarmonics
    for (int nmax = 5; nmax <= maxHarmonics; nmax += 5) {

        std::cout << "\n=== Testing with " << nmax << " harmonics ===\n";

        // Import harmonics for current nmax
        if (importedharmonics < maxHarmonics) {
            freeStokes(importedharmonics);
            importStokesCombined(gravityModels[selectedModel], maxHarmonics);
            importedharmonics = maxHarmonics;
            std::cout << "Harmonics imported.\n";
        }
        
        // Store baseline single-thread Holmes time for speedup calculation
        double baselineTime = 0.0;
        
        // Test Belikov
        {
            std::cout << "Testing Belikov... ";
            std::vector<double> times;
            
            for (int run = 0; run < num_runs; ++run) {
                auto coords = generateRandomCoordinates();
                double radius = coords[0];
                double latitude = coords[1];
                double longitude = coords[2];
                
                std::array<double, 3> result{};
                auto start = std::chrono::high_resolution_clock::now();
                gravityBelikov(radius, latitude, longitude, nmax, result);
                if (isUseImitation)
                    simulate_integrator_and_sofa(27);
                auto end = std::chrono::high_resolution_clock::now();
                double time = std::chrono::duration<double, std::milli>(end - start).count();
                times.push_back(time);
            }
            
            // Calculate statistics
            double avgTime = 0.0;
            for (double t : times) avgTime += t;
            avgTime /= num_runs;
            
            double variance = 0.0;
            for (double t : times) variance += (t - avgTime) * (t - avgTime);
            double stdDev = std::sqrt(variance / num_runs);
            
            double totalTime = avgTime * num_runs;
            
            file << "Belikov,1," << nmax << "," << avgTime << "," << stdDev << "," 
                 << totalTime << ",N/A,N/A\n";
            
            std::cout << "Done. Avg: " << avgTime << " ms\n";
        }
        
        // Test Cunningham
        {
            std::cout << "Testing Cunningham... ";
            std::vector<double> times;
            
            for (int run = 0; run < num_runs; ++run) {
                auto coords = generateRandomCoordinates();
                double radius = coords[0];
                double latitude = coords[1];
                double longitude = coords[2];
                
                std::array<double, 3> result{};
                auto start = std::chrono::high_resolution_clock::now();
                gravityCunningham(radius, latitude, longitude, nmax, result);
                if (isUseImitation)
                    simulate_integrator_and_sofa(27);
                auto end = std::chrono::high_resolution_clock::now();
                double time = std::chrono::duration<double, std::milli>(end - start).count();
                times.push_back(time);
            }
            
            // Calculate statistics
            double avgTime = 0.0;
            for (double t : times) avgTime += t;
            avgTime /= num_runs;
            
            double variance = 0.0;
            for (double t : times) variance += (t - avgTime) * (t - avgTime);
            double stdDev = std::sqrt(variance / num_runs);
            
            double totalTime = avgTime * num_runs;
            
            file << "Cunningham,1," << nmax << "," << avgTime << "," << stdDev << "," 
                 << totalTime << ",N/A,N/A\n";
            
            std::cout << "Done. Avg: " << avgTime << " ms\n";
        }
        
        // Test Holmes with different thread counts (1 to maxThreads)
        for (int threadCount = 1; threadCount <= maxThreads; ++threadCount) {
            std::cout << "Testing Holmes with " << threadCount << " thread(s)... ";
            std::vector<double> times;
            
            using namespace uniorb;
            gravity_stokes GravityStokes(_c, _s, nmax, nmax, EARTH_MU, EARTH_RADIUS);
            GravityStokes.use_concurrency(threadCount);
            
            for (int run = 0; run < num_runs; ++run) {
                auto coords = generateRandomCoordinates();
                double radius = coords[0];
                double latitude = coords[1];
                double longitude = coords[2];
                
                std::array<double, 3> result{};
                auto start = std::chrono::high_resolution_clock::now();
                GravityStokes.get_acceleration(radius, latitude, longitude, result);
                if (isUseImitation)
                    simulate_integrator_and_sofa(27);
                auto end = std::chrono::high_resolution_clock::now();
                double time = std::chrono::duration<double, std::milli>(end - start).count();
                times.push_back(time);
            }
            
            // Calculate statistics
            double avgTime = 0.0;
            for (double t : times) avgTime += t;
            avgTime /= num_runs;
            
            double variance = 0.0;
            for (double t : times) variance += (t - avgTime) * (t - avgTime);
            double stdDev = std::sqrt(variance / num_runs);
            
            double totalTime = avgTime * num_runs;
            
            // Save baseline for speedup calculation
            if (threadCount == 1) {
                baselineTime = avgTime;
            }
            
            // Calculate speedup and efficiency
            double speedup = (threadCount == 1) ? 1.0 : baselineTime / avgTime;
            double efficiency = speedup / threadCount;
            
            std::string algoName = "Holmes" + std::to_string(threadCount);
            file << algoName << "," << threadCount << "," << nmax << "," << avgTime << "," 
                 << stdDev << "," << totalTime << "," << speedup << "," << efficiency << "\n";
            
            std::cout << "Done. Avg: " << avgTime << " ms, Speedup: " << speedup << "\n";
        }
        
        file.flush(); // Ensure data is written periodically
        
        // Progress indicator
        int percentComplete = (int)((double)nmax / maxHarmonics * 100);
        std::cout << "Overall progress: " << percentComplete << "%\n";
    }
    
    file.close();
    std::cout << "\n=== BENCHMARK COMPLETE ===\n";
    std::cout << "Results saved to benchmark_comprehensive.csv\n";
}

// Function to display the benchmarking menu
void displayBenchmarkingMenu() {
    std::cout << "\nBENCHMARKING MODE MENU:\n"
              << "1. Comprehensive Benchmark (All Harmonics & Threads)\n"
              << "2. Back to Main Menu\n"
              << "Enter choice: ";
}

// Main benchmarking mode handler
void runBenchmarkingMode(int& importedharmonics) {
    bool submenu_active = true;
    int num_runs = 100; // Default number of runs per test
    
    while (submenu_active) {
        displayBenchmarkingMenu();
        
        int option;
        std::cin >> option;

        if (std::cin.fail()) {
            std::cin.clear();
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            std::cout << "INVALID INPUT. TRY AGAIN.\n";
            continue;
        }

        switch (option) {
            case 1: {
                std::cout << "Enter number of runs per test (default 100): ";
                std::cin >> num_runs;
                if (std::cin.fail() || num_runs < 1) {
                    std::cin.clear();
                    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                    num_runs = 100;
                    std::cout << "Invalid input. Using default: 100 runs.\n";
                }
                runComprehensiveBenchmark(num_runs, importedharmonics);
                break;
            }
            case 2:
                submenu_active = false;
                break;
            default:
                std::cout << "INVALID OPTION. TRY AGAIN.\n";
                break;
        }
    }
}
