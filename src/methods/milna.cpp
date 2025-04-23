#include "milna.hpp"

std::deque<Eigen::MatrixXd> MilnaSolver::get_acceleration_points(
    int count, double h, Eigen::MatrixXd& P_initial) {
    return acceleration_points(count, h, P_initial);
}

Eigen::MatrixXd MilnaSolver::update_step(
    const Eigen::MatrixXd& P, double h,
    const std::deque<Eigen::MatrixXd>& prev_points) {
    return prev_points[3] +
           ((4 * h / 3.0) *
            (2 * riccati_equation(P) - riccati_equation(prev_points[1]) +
             2 * riccati_equation(prev_points[2])));
}
