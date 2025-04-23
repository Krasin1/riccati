#include "solver.hpp"

#include <Eigen/Dense>
#include <filesystem>
#include <fstream>
#include <iostream>

#include "func.hpp"
extern bool draw;

namespace fs = std::filesystem;

// Наше матричное уравнение Риккати
inline Eigen::MatrixXd RiccatiSolver::riccati_equation(
    const Eigen::MatrixXd& P) {
    return E_ * P * A_ + A_.transpose() * P * E_ + Q_ - E_ * P * BRB_ * P * E_;
}

Result RiccatiSolver::solve(double t0, double t_max, double h,
                            double target_error, int max_steps) {
    Eigen::MatrixXd P = initial_P_;
    Eigen::MatrixXd P_previous = P;
    int step = 0;
    double error = std::numeric_limits<double>::max();

    // если draw == true, то сохраняем точки для графика
    std::vector<double>* points = draw ? new std::vector<double>() : nullptr;

    while (step < max_steps && error > target_error) {
        std::deque<Eigen::MatrixXd> prev_points = acceleration_points(4, h, P);

        P_previous = P;
        for (double t = t0; t < t_max; t += h) {
            progress(target_error, error, t, h, step, t_max, P);

            P = update_step(P, h, prev_points);
            check_nan(P, step, t);

            // prev_points.pop_back();
            // prev_points.push_front(P);
        }
        error = (P - P_previous).norm();
        // значение ошибки на каждом шагу для графика
        if (draw) points->push_back(error);

        step++;
    }

    Result result{P, step, error, points};
    return result;
}

// Найденную матрицу подставляем в уравнение Риккати и выводим в файл
Eigen::MatrixXd RiccatiSolver::verify_solution(const Eigen::MatrixXd& P) {
    Eigen::MatrixXd P_verified = riccati_equation(P);

    fs::path folder_path = "results";
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
std::deque<Eigen::MatrixXd> RiccatiSolver::acceleration_points(
    int count, double h, Eigen::MatrixXd& P_initial) {
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

std::deque<Eigen::MatrixXd> RiccatiSolver::acceleration_points(int count,
                                                               double h) {
    return acceleration_points(count, h, initial_P_);
}
