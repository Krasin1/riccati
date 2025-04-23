#ifndef NYSTROM_HPP
#define NYSTROM_HPP

#include "../core/solver.hpp"

class NystromSolver : public RiccatiSolver {
   public:
    NystromSolver(Eigen::MatrixXd E, Eigen::MatrixXd A, Eigen::MatrixXd B,
                  Eigen::MatrixXd Q, Eigen::MatrixXd initial_P)
        : RiccatiSolver(E, A, B, Q, initial_P) {}

    virtual std::deque<Eigen::MatrixXd> get_acceleration_points(
        int count, double h, Eigen::MatrixXd& P_initial) override;

    virtual Eigen::MatrixXd update_step(
        const Eigen::MatrixXd& P, double h,
        const std::deque<Eigen::MatrixXd>& prev_points) override;
};

#endif
