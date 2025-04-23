#include "runge.hpp"

Eigen::MatrixXd RungeKuttaSolver::update_step(
    const Eigen::MatrixXd& P, double h, const std::deque<Eigen::MatrixXd>&) {
    Eigen::MatrixXd k1 = riccati_equation(P);
    Eigen::MatrixXd k2 = riccati_equation(P + (0.5 * h) * k1);
    Eigen::MatrixXd k3 = riccati_equation(P + (0.5 * h) * k2);
    Eigen::MatrixXd k4 = riccati_equation(P + h * k3);

    return P + (k1 + 2 * k2 + 2 * k3 + k4) * (h / 6.0);
}
