#ifndef FELBERG_H
#define FELBERG_H

#include "../solver.h"

class FelbergSolver : public RiccatiSolver {
   public:
    FelbergSolver(Eigen::MatrixXd E, Eigen::MatrixXd A, Eigen::MatrixXd B,
                     Eigen::MatrixXd Q, Eigen::MatrixXd initial_P)
        : RiccatiSolver(E, A, B, Q, initial_P) {}
    SolverResult solve(double t0, double t_max, double h,
                       double target_error, int max_steps) override;
};

#endif  
