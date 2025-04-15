#include "runge.h"

#include "../../functions/func.h"

extern bool draw;

SolverResult RungeKuttaSolver::solve(double t0, double t_max, double h,
                                     double target_error, int max_steps) {
    Eigen::MatrixXd P = initial_P;
    Eigen::MatrixXd P_previous = P;
    int step = 0;
    double error = std::numeric_limits<double>::max();

    // если draw == true, то сохраняем точки для графика
    std::vector<double>* points = draw ? new std::vector<double>() : nullptr;

    while (step < max_steps && error > target_error) {
        P_previous = P;

        for (double t = t0; t < t_max; t += h) {
            progress(target_error, error, t, h, step, t_max, P);

            Eigen::MatrixXd k1 = riccati_equation(P);
            Eigen::MatrixXd k2 = riccati_equation(P + (0.5 * h) * k1);
            Eigen::MatrixXd k3 = riccati_equation(P + (0.5 * h) * k2);
            Eigen::MatrixXd k4 = riccati_equation(P + h * k3);

            P += (k1 + 2 * k2 + 2 * k3 + k4) * (h / 6.0);

            if (std::isnan(P(0, 0))) {
                std::string message = "Матрица P заполнилась -nan\nШаг : " +
                                      std::to_string(step) +
                                      " | t : " + std::to_string(t) + "\n";
                throw std::runtime_error(message);
            }
        }
        error = (P - P_previous).norm();
        // значение ошибки на каждом шагу для графика
        if (draw) points->push_back(error);

        step++;
    }

    SolverResult result{P, step, error, points};
    return result;
}
