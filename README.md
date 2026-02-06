[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![C++20](https://img.shields.io/badge/C++-20-blue.svg)](https://en.cppreference.com/w/cpp/20)

# 3AKB 🌍 — Three Algorithms for Gravitational Acceleration Benchmarking

**3AKB** is a lightweight, header-only C++20 console application for computing Earth's gravitational acceleration using high-degree spherical harmonic models (EGM96 and EGM2008).

It provides a unified framework to evaluate and benchmark three classical recursive algorithms — Belikov, Cunningham, and Holmes (Stokes) — under identical conditions, including a **multi-threaded parallel implementation** ⚡ of the Holmes algorithm described in our forthcoming paper.

**Key features**
- 🚀 Interactive console menu for quick experiments
- 📍 Single-point evaluation
- 🔄 Random-point and sequential orbit propagation comparison modes
- 📊 Comprehensive multithreaded benchmarking with speedup/efficiency metrics
- 📄 CSV export for post-processing and research workflows
- ✅ No external dependencies

## Table of Contents
- [Quickstart 🚀](#quickstart)
- [Installation 🛠️](#installation)
- [Project structure 📂](#project-structure)
- [Usage 🎮](#usage)
- [Input parameters ⚙️](#input-parameters)
- [Output format 📈](#output-format)
- [Algorithms / API reference 🔧](#algorithms--api-reference)
- [Modes 🧪](#modes)
- [Authors 👥](#authors)
- [Citation 📚](#citation)
- [License 📄](#license)

## Quickstart 🚀

1. Clone the repository:
   ```bash
   git clone https://github.com/timurezy/3AKB.git
   cd 3AKB
   ```
2. Open the Visual Studio solution `3AKB.sln`.
3. Build in **Release** configuration (x64 recommended).
4. Run the executable.

You will be greeted by an interactive menu:
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

**Example**: select EGM2008, set NMAX = 2000, 12 threads → run benchmarking → get detailed CSV with speedup up to ~8.4× ⚡

## Installation 🛠️

### Requirements
- Windows 10/11
- Visual Studio 2019 or 2022 (with C++20 support)
- MSVC toolset

### Build
Open `3AKB.sln` → Build → Build Solution.  
No external libraries or package managers are required.

## Project structure 📂

**Core files**
- `3AKB.cpp` – program entry point and interactive menu
- `alltypes.h` – shared constants, gravity model selection, harmonic storage
- `SingleAlgorithmExecution.h` – single-point evaluation
- `AlgorithmComparison.h` – comparison mode + CSV export
- `AlgorithmBenchmarking.h` – multithreaded benchmarking
- `Simulate.h` – coordinate utilities and random point generation

**Algorithm implementations**
- `Belikov.h` – Belikov method
- `Cunningham.h` – Cunningham method
- `Albert.h` – single- and multi-threaded Holmes (Stokes) implementation

**Data**
- `EGM96.dat`, `EGM2008.dat` – harmonic coefficient files

## Usage 🎮

The program starts with an interactive console menu (see Quickstart).  
Note:
- `NMAX = 0` → central-body (Keplerian) acceleration only
- `NMAX > 0` → harmonic coefficients must be imported first

## Input parameters ⚙️

- **Gravity model**: EGM96 or EGM2008
- **Maximum degree/order (NMAX)**: 0 (Keplerian) to full model resolution
- **Point coordinates**: radius (m), geocentric latitude (°), longitude (°)
- **Number of threads**: 1–maximum supported by hardware (affects parallel Holmes only)
- **Number of runs**: used in comparison and benchmarking modes

## Output format 📈

All results are exported as CSV files in the executable directory.

**Random comparison** – `results_RANDOM.csv`  
Header: `Algorithm,Run,Time (ms),R,Latitude,Longitude,ax,ay,az`  
Includes summary statistics (total time, speedup, efficiency).

**Sequential propagation** – `results_SEQUENTIAL.csv`  
Similar format + final propagated position.

**Benchmarking** – multiple CSV files with detailed timing, speedup, and efficiency metrics.

## Algorithms / API reference 🔧

### Belikov
```cpp
void gravityBelikov(double r, double lat, double lon, int nmax,
                    std::array<double,3>& result);
```

### Cunningham
```cpp
void gravityCunningham(double r, double lat, double lon, int nmax,
                       std::array<double,3>& result);
```

### Holmes / Stokes (parallel-capable) ⚡
```cpp
namespace uniorb {
class gravity_stokes {
public:
    void use_concurrency(int thread_count);
    void get_acceleration(double r, double lat, double lon,
                          int nmax, std::array<double,3>& result);
};
}
```

## Modes 🧪

- **Individual algorithm**: single-point evaluation of selected method
- **Comparison**: random points or simple sequential propagation, with CSV export
- **Benchmarking**: systematic performance tests across degrees and thread counts

## Authors 👥

A. V. Fraerman¹, T. S. Pavlov¹, A. R. Shaykhutdinov¹, P. R. Zapevalin²  
¹ Astro Space Center, Lebedev Physical Institute, Russian Academy of Sciences  
² Sternberg Astronomical Institute, Lomonosov Moscow State University  

Corresponding author: fraerman@asc.rssi.ru

## Citation 📚

If you use this code in your research, please cite our paper (to be updated upon acceptance):

```bibtex
@article{Fraerman2026,
  title = {Multi-threaded parallel implementation of Holmes' algorithm for high-degree spherical harmonic evaluation of Earth's gravitational field},
  author = {Fraerman, A. V. and Pavlov, T. S. and Shaykhutdinov, A. R. and Zapevalin, P. R.},
  journal = {Advances in Space Research},
  year = {2026},
  doi = {...}
}
```

## License 📄

This project is licensed under the MIT License – see the [LICENSE](LICENSE) file for details.
