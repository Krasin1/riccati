#include "felberg.hpp"

Eigen::MatrixXd FelbergSolver::update_step(const Eigen::MatrixXd& P, double h,
                                           const std::deque<Eigen::MatrixXd>&) {
    Eigen::MatrixXd k1 = riccati_equation(P);
    Eigen::MatrixXd k2 = riccati_equation(P + (k1 * (h / 4.0)));
    Eigen::MatrixXd k3 =
        riccati_equation(P + (3 * h / 32.0) * k1 + (9 * h / 32.0) * k2);
    Eigen::MatrixXd k4 =
        riccati_equation(P + (1932 * h / 2197.0) * k1 -
                         (7200 * h / 2197.0) * k2 + (7296 * h / 2197.0) * k3);
    Eigen::MatrixXd k5 =
        riccati_equation(P + (439 * h / 216.0) * k1 - (8 * h) * k2 +
                         (3680 * h / 513.0) * k3 - (845 * h / 4104.0) * k4);

    return P + h * ((25 * k1 / 216.0) + (1408 * k3 / 2565.0) +
                    (2197 * k4 / 4101.0) - (k5 / 5.0));
}
