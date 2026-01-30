# 3AKB

3AKB is a C++ console application for computing Earth gravitational acceleration from spherical harmonic gravity models and evaluating multiple gravity-field algorithms under identical conditions.

The project provides an interactive environment to:

- select a gravity model (EGM96 / EGM2008)
- import harmonic coefficients up to a chosen maximum degree/order
- compute acceleration at a single point
- compare algorithms (random points or sequential propagation)
- benchmark runtime and multithread scaling
- export results to CSV for research workflows

---

Quickstart  
Installation  
Project structure  
Usage  
Input parameters  
Output format  
Algorithms / API reference  
Modes (Individual / Comparison / Benchmarking)  
Authors  
If you use this code in your research, please cite  
License  
Keywords  

C++17, geopotential, spherical harmonics, gravity acceleration, EGM96, EGM2008, multithreading, benchmarking.

---

## 1. Quickstart

1. Open the Visual Studio solution:

```
3AKB.sln
```

2. Build in **Release** mode  
3. Run the executable

You will see an interactive console menu where you can configure parameters and run different modes.

---

## 2. Installation

### Requirements

- Windows
- Visual Studio 2019 / 2022
- C++17 toolset (MSVC)

### Build

Open `3AKB.sln` → Build → Run

No external libraries are required.

---

## 3. Project structure

Core files:

- `3AKB.cpp`  
  Program entry point (`main`) and interactive main menu

- `alltypes.h`  
  Shared constants and global storage:
  - Earth constants (`EARTH_RADIUS`, `EARTH_MU`)
  - gravity model selection
  - harmonic storage arrays
  - import/free helpers

- `SingleAlgorithmExecution.h`  
  Individual algorithm execution mode

- `AlgorithmComparison.h`  
  Algorithm comparison mode + CSV export

- `AlgorithmBenchmarking.h`  
  Benchmarking mode

- `Simulate.h`  
  Coordinate conversions, propagation helpers, random coordinate generation

Algorithm implementations:

- `Belikov.h` → Belikov harmonic summation
- `Cunningham.h` → Cunningham method
- `Albert.h` → Stokes/Holmes implementation (single + multithread)

Data files:

- `EGM96.dat`
- `EGM2008.dat`

---

## 4. Usage

At runtime the program shows:

```
1. RUN INDIVIDUAL ALGORITHM
2. ALGORITHM COMPARISON MODE
3. ALGORITHM BENCHMARKING MODE
4. CHANGE MAX HARMONICS
5. CHANGE INPUT COORDINATES
6. SELECT GRAVITY MODEL
7. IMPORT HARMONICS
8. NUMBER OF THREADS
0. EXIT
```

Notes:

- `NMAX = 0` → Keplerian (central-body) mode
- `NMAX > 0` → spherical harmonics must be imported

---

## 5. Input parameters

### Gravity model
- EGM96
- EGM2008

### Harmonics
- `NMAX` = maximum degree/order
- `NMAX = 0` → Keplerian only

### Coordinates
- radius (meters)
- latitude (degrees)
- longitude (degrees)

### Threads
- affects multithread Stokes solver only
- `threads >= 1`

### Number of runs
Used in comparison and benchmarking modes.

---

## 6. Output format

The program writes CSV files.

### Random comparison

File:

```
results_RANDOM.csv
```

Header:

```
Algorithm,Run,Time (ms),R,Latitude,Longitude,ax,ay,az
```

Includes summary:

- total time per algorithm
- speedup
- efficiency
- final accelerations

---

### Sequential comparison

File:

```
results_SEQUENTIAL.csv
```

Same header + summary:

- total time
- speedup / efficiency
- final propagated coordinates

---

### Benchmarking mode

Creates benchmark CSV files containing:

- algorithm
- threads
- NMAX
- timing statistics
- speedup / efficiency

---

## 7. Algorithms / API reference

### Belikov

```
void gravityBelikov(double r, double lat, double lon, int nmax, std::array<double,3>& result)
```

Computes Cartesian acceleration using Belikov harmonic expansion.

Requires imported harmonics.

---

### Cunningham

```
void gravityCunningham(double r, double lat, double lon, int nmax, std::array<double,3>& result)
```

Uses normalized Legendre recursion and Cunningham summation.

Requires imported harmonics.

---

### Stokes / Holmes (Albert implementation)

Class:

```
uniorb::gravity_stokes
```

Key methods:

```
use_concurrency(thread_count)
get_acceleration(...)
```

Supports automatic multithread splitting of harmonic intervals.

---

## 8. Modes

### Individual algorithm mode

Runs selected solver at current coordinates.

### Comparison mode

Two sub-modes:

1. Random coordinates
2. Sequential propagation

Exports CSV with timings and acceleration vectors.

### Benchmarking mode

Runs repeated experiments for performance studies.

---

## 9. Authors

A. V. Fraerman  
T. S. Pavlov  
A. R. Shaykhutdinov  
P. R. Zapevalin  

Astro Space Center, Lebedev Physical Institute, RAS  
Sternberg Astronomical Institute, Moscow State University  

Corresponding author: fraerman@asc.rssi.ru

---

## 10. If you use this code in your research, please cite:

(Add paper reference here)

---

## 11. License

This project is licensed under the MIT License - see the LICENSE file for details.

---

## 12. Keywords

Earth gravity field  
spherical harmonics  
EGM96  
EGM2008  
orbital mechanics  
multithreading  
benchmarking  
C++17
