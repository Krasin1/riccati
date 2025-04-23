#include "hemming.hpp"

std::deque<Eigen::MatrixXd> HemmingSolver::get_acceleration_points(
    int count, double h, Eigen::MatrixXd& P_initial) {
    return acceleration_points(count, h, P_initial);
}

Eigen::MatrixXd HemmingSolver::update_step(
    const Eigen::MatrixXd& P, double h,
    const std::deque<Eigen::MatrixXd>& prev_points) {
    return (0.5 * (P + prev_points[1])) +
           (h / 48.0) * (119 * riccati_equation(P) -
                         99 * riccati_equation(prev_points[1]) +
                         69 * riccati_equation(prev_points[2]) -
                         17 * riccati_equation(prev_points[3]));
}
