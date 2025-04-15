#include "adams.h"

#include "../../functions/func.h"

extern bool draw;

SolverResult AdamsSolver::solve(double t0, double t_max, double h,
                                double target_error, int max_steps) {
    Eigen::MatrixXd P = initial_P;
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

            P = P + (h / 24.0) * ((55 * riccati_equation(P)) -
                                  (59 * riccati_equation(prev_points[1])) +
                                  (37 * riccati_equation(prev_points[2])) -
                                  (9 * riccati_equation(prev_points[3])));

            if (std::isnan(P(0, 0))) {
                std::string message = "Матрица P заполнилась -nan\nШаг : " +
                                      std::to_string(step) +
                                      " | t : " + std::to_string(t) + "\n";
                throw std::runtime_error(message);
            }

            prev_points.pop_back();
            prev_points.push_front(P);
        }
        error = (P - P_previous).norm();
        // значение ошибки на каждом шагу для графика
        if (draw) points->push_back(error);

        step++;
    }

    SolverResult result{P, step, error, points};
    return result;
}
