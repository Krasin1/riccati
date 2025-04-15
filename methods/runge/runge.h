#ifndef RUNGE_KUTTA_H
#define RUNGE_KUTTA_H

#include "../solver.h"

class RungeKuttaSolver : public RiccatiSolver {
   public:
    RungeKuttaSolver(Eigen::MatrixXd E, Eigen::MatrixXd A, Eigen::MatrixXd B,
                     Eigen::MatrixXd Q, Eigen::MatrixXd initial_P)
        : RiccatiSolver(E, A, B, Q, initial_P) {}
    SolverResult solve(double t0, double t_max, double h,
                       double target_error, int max_steps) override;
};

#endif  // RUNGE_KUTTA_H
