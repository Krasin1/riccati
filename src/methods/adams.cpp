#include "adams.hpp"

std::deque<Eigen::MatrixXd> AdamsSolver::get_acceleration_points(
    int count, double h, Eigen::MatrixXd& P_initial) {
    return acceleration_points(count, h, P_initial);
}

Eigen::MatrixXd AdamsSolver::update_step(
    const Eigen::MatrixXd& P, double h,
    const std::deque<Eigen::MatrixXd>& prev_points) {
    return P + (h / 24.0) * ((55 * riccati_equation(P)) -
                             (59 * riccati_equation(prev_points[1])) +
                             (37 * riccati_equation(prev_points[2])) -
                             (9 * riccati_equation(prev_points[3])));
}
