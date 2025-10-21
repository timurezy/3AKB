#pragma once
#include <cmath>
#include <vector>
#include <random>
#include <chrono>
#include <iostream>
#include <array>

// Имитационная функция с настраиваемым количеством итераций
void simulate_integrator_and_sofa(int num_iterations) {
    const int matrix_size = 3; // Размер матриц 3x3 для преобразований СК
    const int array_size = 10000; // Размер массива для имитации EOP-данных

    // Генерация случайных данных для предотвращения оптимизаций
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(-1.0, 1.0);

    // Массив для имитации EOP-данных
    std::vector<double> eop_data(array_size);
    for (auto& val : eop_data) {
        val = dist(gen);
    }

    // Матрицы для преобразований СК
    std::vector<std::vector<double>> matrix_a(matrix_size, std::vector<double>(matrix_size));
    std::vector<std::vector<double>> matrix_b(matrix_size, std::vector<double>(matrix_size));
    std::vector<std::vector<double>> matrix_result(matrix_size, std::vector<double>(matrix_size));

    // Заполнение матриц случайными значениями
    for (int i = 0; i < matrix_size; ++i) {
        for (int j = 0; j < matrix_size; ++j) {
            matrix_a[i][j] = dist(gen);
            matrix_b[i][j] = dist(gen);
        }
    }

    // Имитация вычислений
    volatile double result = 0.0; // Предотвращает оптимизацию
    for (int iter = 0; iter < num_iterations; ++iter) {
        // 1. Тригонометрические вычисления (имитация SOFA)
        for (int i = 0; i < array_size; ++i) {
            result += std::sin(eop_data[i]) * std::cos(eop_data[i]);
        }

        // 2. Матричные операции (имитация преобразований СК)
        for (int i = 0; i < matrix_size; ++i) {
            for (int j = 0; j < matrix_size; ++j) {
                matrix_result[i][j] = 0.0;
                for (int k = 0; k < matrix_size; ++k) {
                    matrix_result[i][j] += matrix_a[i][k] * matrix_b[k][j];
                }
            }
        }

        // 3. Имитация чтения EOP-данных (доступ к памяти)
        for (int i = 0; i < array_size; ++i) {
            eop_data[i] = dist(gen) * result; // Зависимость данных
        }
    }

    // Используем результат, чтобы компилятор не удалил вычисления
    if (result > 0) {
        matrix_a[0][0] = result; // Простое использование результата
    }
}




inline std::array<double, 3> generateRandomCoordinates() {
    static std::random_device rd;
    static std::mt19937 gen(rd());

    std::uniform_real_distribution<double> distRadius(6000000.0, 7000000.0);
    std::uniform_real_distribution<double> distLat(-90.0, 90.0);
    std::uniform_real_distribution<double> distLon(-180.0, 180.0);

    double r = distRadius(gen);
    double lat = distLat(gen);
    double lon = distLon(gen);

    return { r, lat, lon };
}
