#include "inglend.hpp"

Eigen::MatrixXd InglendSolver::update_step(const Eigen::MatrixXd& P, double h,
                                           const std::deque<Eigen::MatrixXd>&) {
    Eigen::MatrixXd k1 = riccati_equation(P);
    Eigen::MatrixXd k2 = riccati_equation(P + (k1 * (h / 2.0)));
    Eigen::MatrixXd k3 = riccati_equation(P + (h / 4.0) * (k1 + k2));
    Eigen::MatrixXd k4 = riccati_equation(P + h * (2 * k3 - k2));

    return P + (h / 6.0) * (k1 + 4 * k3 + k4);
}
