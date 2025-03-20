#ifndef HEMMING_H
#define HEMMING_H 

#include "../solver.h"

class HemmingSolver : public RiccatiSolver {
   public:
    HemmingSolver(Eigen::MatrixXd E, Eigen::MatrixXd A, Eigen::MatrixXd B,
                Eigen::MatrixXd Q, Eigen::MatrixXd initial_P)
        : RiccatiSolver(E, A, B, Q, initial_P) {}
    SolverResult solve(double t0, double t_max, double h,
                       double target_error) override;
};

#endif
