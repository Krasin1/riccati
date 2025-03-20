#include "solver.h"

#include <Eigen/Dense>
#include <filesystem>
#include <fstream>
#include <iostream>

namespace fs = std::filesystem;

// Наше матричное уравнение Риккати
Eigen::MatrixXd RiccatiSolver::riccati_equation(const Eigen::MatrixXd& P) {
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
std::deque<Eigen::MatrixXd> RiccatiSolver::acceleration_points(int count,
                                                               double h) {
    std::deque<Eigen::MatrixXd> result;
    Eigen::MatrixXd k1, k2, k3, k4;

    result.push_front(initial_P);
    for (int i = 1; i < count; i++) {
        k1 = riccati_equation(initial_P);
        k2 = riccati_equation(initial_P + (k1 * h / 2.0));
        k3 = riccati_equation(initial_P + (k2 * h / 2.0));
        k4 = riccati_equation(initial_P + (k3 * h));

        initial_P += (k1 + 2.0 * k2 + 2.0 * k3 + k4) * (h / 6.0);
        result.push_front(initial_P);
    }

    return result;
}
