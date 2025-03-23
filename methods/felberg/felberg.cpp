#include "felberg.h"

#include "../../functions/func.h"

extern bool draw;

SolverResult FelbergSolver::solve(double t0, double t_max, double h,
                                  double target_error) {
    Eigen::MatrixXd P = initial_P;
    Eigen::MatrixXd P_previous = P;
    int step = 0;
    double error = std::numeric_limits<double>::max();

    // если draw == true, то сохраняем точки для графика
    std::vector<double>* points = draw ? new std::vector<double>() : nullptr;

    while (step < 200 && error > target_error) {
        P_previous = P;

        for (double t = t0; t < t_max; t += h) {
            progress(target_error, error, t, h, step, t_max, P);

            Eigen::MatrixXd k1 = riccati_equation(P);
            Eigen::MatrixXd k2 = riccati_equation(P + (k1 * (h / 4.0)));
            Eigen::MatrixXd k3 =
                riccati_equation(P + (3 * h / 32.0) * k1 + (9 * h / 32.0) * k2);
            Eigen::MatrixXd k4 = riccati_equation(P + (1932 * h / 2197.0) * k1 -
                                                  (7200 * h / 2197.0) * k2 +
                                                  (7296 * h / 2197.0) * k3);
            Eigen::MatrixXd k5 = riccati_equation(
                P + (439 * h / 216.0) * k1 - (8 * h) * k2 +
                (3680 * h / 513.0) * k3 - (845 * h / 4104.0) * k4);

            P += h * ((25 * k1 / 216.0) + (1408 * k3 / 2565.0) +
                      (2197 * k4 / 4101.0) - (k5 / 5.0));

            if (std::isnan(P(0, 0))) {
                std::string message = "Матрица P заполнилась -nan\nШаг : " +
                                      std::to_string(step) +
                                      " | t : " + std::to_string(t) + "\n";
                throw std::runtime_error(message);
            }
        }

        // значение ошибки на каждом шагу для графика
        if (draw) points->push_back(error);

        error = (P - P_previous).norm();
        step++;
    }

    SolverResult result{P, step, error, points};
    return result;
}
