#include "milna.h"

#include "../../functions/func.h"

extern bool draw;
extern bool manual;

SolverResult MilnaSolver::solve(double t0, double t_max, double h,
                                double target_error) {
    Eigen::MatrixXd P = initial_P;
    Eigen::MatrixXd P_previous = P;

    int step = 0;
    double error = std::numeric_limits<double>::max();

    // если draw == true, то сохраняем точки для графика
    std::vector<double>* points = draw ? new std::vector<double>() : nullptr;

    // получаем 4 разгонные точки
    std::deque<Eigen::MatrixXd> prev_points = acceleration_points(4, h);

    while (step < 200 && error > target_error) {
        P_previous = P;

        for (double t = t0; t < t_max; t += h) {
            progress(target_error, error, t, h, step, t_max, P);

            P = prev_points[3] +
                ((4 * h / 3.0) * (2 * riccati_equation(prev_points[0]) -
                                  riccati_equation(prev_points[1]) +
                                  2 * riccati_equation(prev_points[2])));

            if (std::isnan(P(0, 0))) {
                std::string message = "Матрица P заполнилась -nan\nШаг : " +
                                      std::to_string(step) +
                                      " | t : " + std::to_string(t) + "\n";
                throw std::runtime_error(message);
            }

            prev_points.pop_back();
            prev_points.push_front(P);
        }
        // значение ошибки на каждом шагу для графика
        if (draw) points->push_back(error);

        error = (P - P_previous).norm();
        step++;
    }

    SolverResult result{P, step, error, points};
    return result;
}
