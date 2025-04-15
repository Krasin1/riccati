#include "solver.h"

#include <Eigen/Dense>
#include <filesystem>
#include <fstream>
#include <iostream>

namespace fs = std::filesystem;

// Наше матричное уравнение Риккати
inline Eigen::MatrixXd RiccatiSolver::riccati_equation(const Eigen::MatrixXd& P) {
    return E * P * A + A.transpose() * P * E + Q - E * P * BRB * P * E;
}

// Найденную матрицу подставляем в уравнение Риккати и выводим в файл
Eigen::MatrixXd RiccatiSolver::verify_solution(const Eigen::MatrixXd& P) {
    Eigen::MatrixXd P_verified = riccati_equation(P);

    fs::path folder_path = "Results";
    fs::path P_verified_file_path = folder_path / "verified.txt";

    if (!fs::exists(folder_path)) fs::create_directory(folder_path);
    std::ofstream verified(P_verified_file_path);

    if (verified.is_open()) {
        verified << std::fixed << std::setprecision(7);
        verified << P_verified << '\n';
        verified.close();
    } else {
        std::cout << "Unable to open file\n";
    }
    return riccati_equation(P);
}

// С помощью метода Рунге-Кутты получаем начальные точки для некоторых методов
std::deque<Eigen::MatrixXd> RiccatiSolver::acceleration_points(int count, double h, Eigen::MatrixXd& P_initial) {
    std::deque<Eigen::MatrixXd> result;
    Eigen::MatrixXd k1, k2, k3, k4;

    for (int i = 0; i < count; i++) {
        k1 = h * riccati_equation(P_initial);
        k2 = h * riccati_equation(P_initial + (k1 / 2.0));
        k3 = h * riccati_equation(P_initial + (k2 / 2.0));
        k4 = h * riccati_equation(P_initial + (k3));

        P_initial += ((k1 + 2.0 * k2 + 2.0 * k3 + k4) / 36.0);
        result.push_front(P_initial);
    }

    return result;
}

std::deque<Eigen::MatrixXd> RiccatiSolver::acceleration_points(int count, double h) {
    return acceleration_points(count, h, initial_P);
}

