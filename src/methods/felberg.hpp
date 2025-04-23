#ifndef FELBERG_HPP
#define FELBERG_HPP

#include "../core/solver.hpp"

class FelbergSolver : public RiccatiSolver {
   public:
    FelbergSolver(Eigen::MatrixXd E, Eigen::MatrixXd A, Eigen::MatrixXd B,
                  Eigen::MatrixXd Q, Eigen::MatrixXd initial_P)
        : RiccatiSolver(E, A, B, Q, initial_P) {}

    Eigen::MatrixXd update_step(
        const Eigen::MatrixXd& P, double h,
        const std::deque<Eigen::MatrixXd>& prev_points) override;
};

#endif
