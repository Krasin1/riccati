#include "nystrom.hpp"

std::deque<Eigen::MatrixXd> NystromSolver::get_acceleration_points(
    int count, double h, Eigen::MatrixXd& P_initial) {
    return acceleration_points(count, h, P_initial);
}

Eigen::MatrixXd NystromSolver::update_step(
    const Eigen::MatrixXd& P, double h,
    const std::deque<Eigen::MatrixXd>& prev_points) {
    return prev_points[1] + (h / 3.0) * (8 * riccati_equation(P) -
                                         5 * riccati_equation(prev_points[1]) +
                                         4 * riccati_equation(prev_points[2]) -
                                         riccati_equation(prev_points[3]));
}
