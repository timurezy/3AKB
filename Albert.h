#pragma once
#include "alltypes.h"
#include <vector>
#include <array>
#include <cmath>
#include <functional>
#include <fstream>
#include <sstream>
#include <iostream>
#include <future>
#include <filesystem>

namespace uniorb {

    class gravity_stokes {
    public:
        struct Spherical {
            double Vlon = 0;
            double Vtet = 0;
            double Vr = 0;
        };

        static constexpr double PI1_2 = 1.5707963267948966;

    public:
        gravity_stokes() :
            _c(),
            _s(),
            _order(0),
            _degree(0),
            _gm(0),
            _radius(0),
            _intervals()
        {
            use_concurrency(1);
        }

        gravity_stokes(std::vector<std::vector<double>>&& C, std::vector<std::vector<double>>&& S, size_t Order, size_t Degree, double GM, double Radius) :
            _c(C),
            _s(S),
            _order(Order),
            _degree(Degree),
            _gm(GM),
            _radius(Radius),
            _intervals()
        {
            use_concurrency(1);
        }

        gravity_stokes(const gravity_stokes& GravityStokes) = default;
        gravity_stokes& operator=(const gravity_stokes& GravityStokes) = default;
        gravity_stokes(gravity_stokes&& GravityStokes) = default;
        gravity_stokes& operator=(gravity_stokes&& GravityStokes) = default;
        virtual ~gravity_stokes() = default;

        void use_concurrency(size_t ThreadCount) {
            ThreadCount = std::max((size_t)1, ThreadCount);
            if (ThreadCount == 1) {
                _method = [](const gravity_stokes& GS, double R, double Lat, double Long, std::array<double, 3>& A) { GS.get_acceleration_one_thread(R, Lat, Long, A); };
            }
            else {
                compute_intervals(ThreadCount);
                _method = [](const gravity_stokes& GS, double R, double Lat, double Long, std::array<double, 3>& A) { GS.get_acceleration_multithread(R, Lat, Long, A); };
            }
        }

        void get_acceleration(double R, double Lat, double Long, std::array<double, 3>& A) const {
            _method(*this, R, Lat, Long, A);
        }

        /*//void import(const std::filesystem::path& FilePath, size_t Order, size_t Degree) 
        void import(const std::string& FilePath, size_t Order, size_t Degree) {
            // Prepare to read data.
            std::cout << "INFO  | Start reading gravity_stokes data from file " << FilePath << std::endl;
            std::ifstream file(FilePath);
            std::string line;
            size_t max_order;
            size_t max_degree;

            if (file.good()) {
                // Read GM
                if (std::getline(file, line)) {
                    _gm = std::stod(line);
                }
                else {
                    std::cout << "ERROR | Failed to read GM, " << FilePath << std::endl;
                    file.close();
                }

                // Read equatorial radius
                if (std::getline(file, line)) {
                    _radius = std::stod(line);
                }
                else {
                    std::cout << "ERROR | Failed to read equatorial radius, " << FilePath << std::endl;
                    file.close();
                }

                // Read the order
                if (std::getline(file, line)) {
                    max_order = std::stoi(line);
                }
                else {
                    std::cout << "ERROR | Failed to read max order, " << FilePath << std::endl;
                    file.close();
                    return;
                }

                // Read the degree
                if (std::getline(file, line)) {
                    max_degree = std::stoi(line);
                }
                else {
                    std::cout << "ERROR | Failed to read max degree, " << FilePath << std::endl;
                    file.close();
                    return;
                }

                // Skip the line
                std::getline(file, line);

                // Check order and degree
                if (max_order < Order) {
                    std::cout << "WARN  | Required order is less than the given one. " << max_order << " will be used, " << FilePath << std::endl;
                    _order = max_order;
                }
                if (max_degree < Degree) {
                    std::cout << "WARN  | Required degree is less than the given one. " << max_degree << " will be used, " << FilePath << std::endl;
                    Degree = max_degree;
                }
                if (Order < Degree) {
                    std::cout << "WARN  | Order can not be less than degree. " << Order << " will be used as degree value, " << FilePath << std::endl;
                    Order = Degree;
                }

                // Initialize vectors for stokes
                _c = std::vector<std::vector<double>>(Order + 1);
                _s = std::vector<std::vector<double>>(Order + 1);
                auto IsLoaded = std::vector<std::vector<bool>>(Order + 1);
                size_t degree = Degree;
                size_t total = 0;
                for (size_t n = 0; n <= Order; ++n) {
                    degree = std::min(n, Degree);
                    _c[n] = std::vector<double>(degree + 1);
                    _s[n] = std::vector<double>(degree + 1);
                    IsLoaded[n] = std::vector<bool>(degree + 1);
                    total += (n + 1);
                }

                // Read the stokes
                size_t n, m;
                double c, s;
                while (total > 0 && std::getline(file, line)) {
                    std::istringstream iss(line);
                    if (iss >> n >> m >> c >> s) {
                        if (n <= Order && m <= Degree && IsLoaded[n][m] == false) {
                            _c[n][m] = c;
                            _s[n][m] = s;
                            IsLoaded[n][m] = true;
                            --total;
                        }
                    }
                    else {
                        std::cout << "ERROR | Failed to read C and S data, " << FilePath << std::endl;
                        file.close();
                        return;
                    }
                }
                // If total > 0 - not all data are present in file
                if (total > 0) {
                    std::cout << "WARN  | C and S are missing: " << total << " values. 0 will be used instead of them," << FilePath << std::endl;
                }

                // Return gravity_stokes
                _order = Order;
                _degree = Degree;

                std::cout << "INFO  | Graivty field loaded successfully from file " << FilePath << std::endl;
                file.close();
            }
            else {
                std::cout << "ERROR | Failed to load gravity data from file " << FilePath << std::endl;
            }
        }
*/



void import(const std::string& FilePath, size_t Order, size_t Degree) {
    std::cout << "INFO  | Start reading gravity_stokes data from file " << FilePath << std::endl;

    _gm = EARTH_MU;
    _radius = EARTH_RADIUS;
    _order = Order;
    _degree = Degree;

    // Инициализация векторов для коэффициентов
    _c = std::vector<std::vector<double>>(Order + 1);
    _s = std::vector<std::vector<double>>(Order + 1);
    for (size_t n = 0; n <= Order; ++n) {
        _c[n] = std::vector<double>(Degree + 1, 0.0);
        _s[n] = std::vector<double>(Degree + 1, 0.0);
    }

    std::ifstream file(FilePath);
    if (!file.is_open()) {
        std::cout << "ERROR | Failed to open file " << FilePath << std::endl;
        return;
    }

    size_t n, m;
    double c, s;
    std::string line;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        if (iss >> n >> m >> c >> s) {
            if (n <= Order && m <= Degree) {
                _c[n][m] = c;
                _s[n][m] = s;
            
        }
    }

    std::cout << "INFO  | Gravity field loaded successfully from file " << FilePath << std::endl;
}




    protected:
        std::vector<std::vector<double>> _c;
        std::vector<std::vector<double>> _s;
        size_t _degree;
        size_t _order;
        double _gm;
        double _radius;
        std::function<void(const gravity_stokes&, double, double, double, std::array<double, 3>&)> _method;
        std::vector< std::pair<size_t, size_t> > _intervals;

    private:
        void get_acceleration_one_thread(double R, double Lat, double Long, std::array<double, 3>& A) const {
            // Support values
            double t = std::cos(PI1_2 - Lat);
            double u = std::sin(PI1_2 - Lat);

            double coslon = std::cos(Long);
            double sinlon = std::sin(Long);

            // Check R
            if (R <= 0) {
                std::cout << "gravity_stokes | _get_acceleration | R must be a positive value." << std::endl;
                return;
            }
            double q0 = _radius / R;

            // Pre-caclulation of Pmm
            std::vector<double> Q = std::vector<double>(_order + 1);
            std::vector<double> um = std::vector<double>(_degree + 1);
            Q[0] = 1.;
            um[0] = 1.;
            if (_order > 0) Q[1] = 1.73205080756887730;
            if (_degree > 0) um[1] = u;
            for (int i = 2; i <= _order; i++) {
                Q[i] = Q[i - 1] * std::sqrt((2. * i + 1.) / (2. * i));
            }
            for (int i = 2; i <= _degree; i++) {
                um[i] = um[i - 1] * u; // Here is the modified method
            }

            // Pre-caclulation of q
            std::vector<double> q = std::vector<double>(_order + 1);
            q[0] = 1.;
            for (int i = 1; i <= _order; i++)
                q[i] = q[i - 1] * q0;

            // Initialize XmC and XmS coefficients and cos/sin lambda	
            double Vlon = 0.;
            double Vtet = 0.;
            double Vr = 0.;

            // Compute
            auto Spherical = compute_partial(q, Q, um, Long, coslon, sinlon, t, u, 0, _degree);

            // Final calulcation of V values
            double GMr = _gm / R;
            Vlon = -GMr * Spherical.Vlon;
            Vtet = GMr * Spherical.Vtet;
            Vr = -GMr / R * Spherical.Vr;

            // Calculation of acceleration
            A[0] = u * coslon * Vr - t / R * coslon * Vtet - sinlon / u / R * Vlon;
            A[1] = u * sinlon * Vr - t / R * sinlon * Vtet + coslon / u / R * Vlon;
            A[2] = t * Vr + u / R * Vtet;
        }

        void get_acceleration_multithread(double R, double Lat, double Long, std::array<double, 3>& A) const {
            // Support values
            double t = std::cos(PI1_2 - Lat);
            double u = std::sin(PI1_2 - Lat);

            double coslon = std::cos(Long);
            double sinlon = std::sin(Long);

            // Check R
            if (R <= 0) {
                std::cout << "gravity_stokes | _get_acceleration_multithread | R must be a positive value." << std::endl;
                return;
            }
            double q0 = _radius / R;

            // Pre-caclulation of Pmm
            std::vector<double> Q = std::vector<double>(_order + 1);
            std::vector<double> um = std::vector<double>(_degree + 1);
            Q[0] = 1.;
            um[0] = 1.;
            if (_order > 0) Q[1] = 1.73205080756887730;
            if (_degree > 0) um[1] = u;
            for (int i = 2; i <= _order; i++) {
                Q[i] = Q[i - 1] * std::sqrt((2. * i + 1.) / (2. * i));
            }
            for (int i = 2; i <= _degree; i++) {
                um[i] = um[i - 1] * u; // Here is the modified method
            }

            // Pre-caclulation of q
            std::vector<double> q = std::vector<double>(_order + 1);
            q[0] = 1.;
            for (int i = 1; i <= _order; i++) {
                q[i] = q[i - 1] * q0;
            }

            // Parallel computing - compute intervals.
            auto Partial = std::vector< std::future<Spherical> >();
            Partial.reserve(_intervals.size());

            // Launch tasks in a thread pool.
            auto f = [&, this](size_t M1, size_t M2) -> Spherical {
                return this->compute_partial(q, Q, um, Long, coslon, sinlon, t, u, M1, M2);
                };

            for (size_t i = 0; i < _intervals.size(); ++i) {
                Partial.emplace_back(std::async(std::launch::async, f, _intervals[i].first, _intervals[i].second));
            }

            // Collect result.
            auto Sum = Spherical();
            auto Part = Spherical();
            for (size_t i = 0; i < _intervals.size(); ++i) {
                Part = Partial[i].get();
                Sum.Vlon += Part.Vlon;
                Sum.Vtet += Part.Vtet;
                Sum.Vr += Part.Vr;
            }

            // Final calulcation of V values
            double GMr = _gm / R;
            Sum.Vlon = -GMr * Sum.Vlon;
            Sum.Vtet = GMr * Sum.Vtet;
            Sum.Vr = -GMr / R * Sum.Vr;

            // Calculation of acceleration
            A[0] = u * coslon * Sum.Vr - t / R * coslon * Sum.Vtet - sinlon / u / R * Sum.Vlon;
            A[1] = u * sinlon * Sum.Vr - t / R * sinlon * Sum.Vtet + coslon / u / R * Sum.Vlon;
            A[2] = t * Sum.Vr + u / R * Sum.Vtet;
        }

        Spherical compute_partial(const std::vector<double>& q, std::vector<double>& Q, std::vector<double>& um, double Long, double coslon, double sinlon, double t, double u, size_t M1, size_t M2) const {
            double XmC, XmS, XmCtet, XmStet, XmCr, XmSr, Qnm, dQnm, Qn_1m, Qn_2m, anm, bnm, qc, qs, tmp;
            double sinl = std::sin(M1 * Long);
            double cosl = std::cos(M1 * Long);
            double sinl_old = 0;
            double cosl_old = 1;
            Spherical Result;

            // Calculation of potential // Outer loop
            for (size_t m = M1; m <= M2; m++) {
                // Get ready n = m loop
                Qnm = Q[m];
                dQnm = -t * m * Qnm / u;
                Qn_1m = Qnm;
                Qn_2m = 0.;
                qc = q[m] * _c[m][m];
                qs = q[m] * _s[m][m];

                tmp = Qnm * um[m];
                XmC = qc * tmp;
                XmS = qs * tmp;
                tmp = dQnm * um[m];
                XmCtet = qc * tmp;
                XmStet = qs * tmp;
                XmCr = (m + 1.) * XmC;
                XmSr = (m + 1.) * XmS;

                // Inner loop // m + 1 beacuse m is included already
                for (size_t n = m + 1; n <= _order; n++) {
                    tmp = (2. * n + 1.) / (n - m) / (n + m);
                    anm = std::sqrt((2. * n - 1.) * tmp);
                    bnm = std::sqrt(tmp * (n + m - 1.) * (n - m - 1.) / (2. * n - 3.));

                    Qnm = anm * t * Qn_1m - bnm * Qn_2m;
                    dQnm = -(n * t * Qnm - (2. * n + 1) * Qn_1m / anm) / u;

                    qc = q[n] * _c[n][m];
                    qs = q[n] * _s[n][m];

                    tmp = Qnm * um[m];
                    XmC += qc * tmp;
                    XmS += qs * tmp;
                    XmCr += (n + 1) * qc * tmp;
                    XmSr += (n + 1) * qs * tmp;
                    tmp = dQnm * um[m];
                    XmCtet += qc * tmp;
                    XmStet += qs * tmp;

                    Qn_2m = Qn_1m;
                    Qn_1m = Qnm;
                }

                // Add new values to sums in V
                Result.Vlon += m * (sinl * XmC - cosl * XmS);
                Result.Vtet += cosl * XmCtet + sinl * XmStet;
                Result.Vr += cosl * XmCr + sinl * XmSr;
                // Prepare for the next loop
                cosl_old = cosl * coslon - sinl * sinlon;
                sinl_old = sinl * coslon + cosl * sinlon;

                // Save cos and sin for recursion
                cosl = cosl_old;
                sinl = sinl_old;
            }
            return Result;
        }

        void compute_intervals(size_t Count) {
            _intervals.clear();

            // Compute total.
            size_t Total = (_degree + 2 * (_order + 1 + _order * _degree) - _degree * _degree) / 2;

            // Find average values in interval.
            size_t Average = std::max((Total / Count), (size_t)1);

            // Find intervals.
            size_t M1 = 0;
            size_t M2 = 0;
            size_t Sum = 0;
            while (_intervals.size() != Count) {
                Sum += _order - M2;
                if (Sum > Average) {
                    _intervals.emplace_back(std::pair<size_t, size_t>(M1, M2));
                    M1 = M2 + 1;
                    Sum = 0;
                }
                if (M2 == _degree) {
                    _intervals.emplace_back(std::pair<size_t, size_t>(M1, M2));
                    break;
                }
                ++M2;
            }
        }

    };

};