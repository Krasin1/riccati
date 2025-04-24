#include "include/solver.hpp"

#include <Eigen/Dense>
#include <deque>
#include <filesystem>
#include <fstream>
#include <iostream>

#include "include/utils.hpp"

// Наше матричное уравнение Риккати
inline Eigen::MatrixXd RiccatiSolver::riccati_equation(
    const Eigen::MatrixXd& P) {
    return E_ * P * A_ + A_.transpose() * P * E_ + Q_ - E_ * P * BRB_ * P * E_;
}

Result RiccatiSolver::solve(Config cfg) {
    Eigen::MatrixXd P = initial_P_;
    Eigen::MatrixXd P_previous = P;
    int step = 0;
    double error = std::numeric_limits<double>::max();

    // если draw == true, то сохраняем точки для графика
    std::vector<double>* points =
        cfg.draw ? new std::vector<double>() : nullptr;

    std::chrono::time_point<std::chrono::system_clock> begin, end;

    while (step < cfg.max_steps && error > cfg.target_error) {
        std::deque<Eigen::MatrixXd> prev_points =
            acceleration_points(4, cfg.h, P);

        if (cfg.step_time) begin = std::chrono::system_clock::now();

        P_previous = P;
        for (double t = cfg.t0; t < cfg.t_max; t += cfg.h) {
            progress(cfg, error, step, t, P);

            P = update_step(P, cfg.h, prev_points);
            check_nan(P, step, t);

            // prev_points.pop_back();
            // prev_points.push_front(P);
        }
        error = (P - P_previous).norm();

        // значение ошибки на каждом шагу для графика
        if (cfg.draw) points->push_back(error);

        if (cfg.step_time) {
            end = std::chrono::system_clock::now();
            auto duration =
                std::chrono::duration_cast<std::chrono::milliseconds>(end -
                                                                      begin);
            std::cout << duration.count() << "ms\n";
        }

        step++;
    }

    Result result{P, step, error, points};
    return result;
}

// Найденную матрицу подставляем в уравнение Риккати и выводим в файл
Eigen::MatrixXd RiccatiSolver::verify_solution(const Eigen::MatrixXd& P) {
    namespace fs = std::filesystem;
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

// С помощью метода Рунге-Кутты получаем начальные точки для некоторых
// методов
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
