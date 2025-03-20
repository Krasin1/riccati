#ifndef MILNA_H
#define MILNA_H

#include "../solver.h"

class MilnaSolver : public RiccatiSolver {
   public:
    MilnaSolver(Eigen::MatrixXd E, Eigen::MatrixXd A, Eigen::MatrixXd B,
                Eigen::MatrixXd Q, Eigen::MatrixXd initial_P)
        : RiccatiSolver(E, A, B, Q, initial_P) {}
    SolverResult solve(double t0, double t_max, double h,
                       double target_error) override;
};

#endif
