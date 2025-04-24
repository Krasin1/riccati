#ifndef METHODS_HPP
#define METHODS_HPP

#include "solver.hpp"

template <typename Method>
class Solver : public RiccatiSolver {
   public:
    Solver(Eigen::MatrixXd& E, Eigen::MatrixXd& A, Eigen::MatrixXd& B,
           Eigen::MatrixXd& Q, Eigen::MatrixXd& initial_P)
        : RiccatiSolver(E, A, B, Q, initial_P) {}

   private:
    Eigen::MatrixXd update_step(
        const Eigen::MatrixXd& P, double h,
        const std::deque<Eigen::MatrixXd>& prev) override {
        return Method::step(*this, P, h, prev);
    };
};

struct RungeKutta {
    static Eigen::MatrixXd step(RiccatiSolver& solver, const Eigen::MatrixXd& P,
                                double h, const std::deque<Eigen::MatrixXd>&) {
        Eigen::MatrixXd k1, k2, k3, k4;
        k1.noalias() = solver.riccati_equation(P);
        k2.noalias() = solver.riccati_equation(P + (0.5 * h) * k1);
        k3.noalias() = solver.riccati_equation(P + (0.5 * h) * k2);
        k4.noalias() = solver.riccati_equation(P + h * k3);

        return P + (k1 + 2 * k2 + 2 * k3 + k4) * (h / 6.0);
    }
};

struct Adams {
    static Eigen::MatrixXd step(RiccatiSolver& solver, const Eigen::MatrixXd& P,
                                double h,
                                const std::deque<Eigen::MatrixXd>& prev) {
        return P + (h / 24.0) * ((55 * solver.riccati_equation(P)) -
                                 (59 * solver.riccati_equation(prev[1])) +
                                 (37 * solver.riccati_equation(prev[2])) -
                                 (9 * solver.riccati_equation(prev[3])));
    }
};

struct Milna {
    static Eigen::MatrixXd step(RiccatiSolver& solver, const Eigen::MatrixXd& P,
                                double h,
                                const std::deque<Eigen::MatrixXd>& prev) {
        return prev[3] +
               ((4 * h / 3.0) * (2 * solver.riccati_equation(P) -
                                 solver.riccati_equation(prev[1]) +
                                 2 * solver.riccati_equation(prev[2])));
    }
};

struct Nystrom {
    static Eigen::MatrixXd step(RiccatiSolver& solver, const Eigen::MatrixXd& P,
                                double h,
                                const std::deque<Eigen::MatrixXd>& prev) {
        return prev[1] + (h / 3.0) * (8 * solver.riccati_equation(P) -
                                      5 * solver.riccati_equation(prev[1]) +
                                      4 * solver.riccati_equation(prev[2]) -
                                      solver.riccati_equation(prev[3]));
    }
};

struct Hemming {
    static Eigen::MatrixXd step(RiccatiSolver& solver, const Eigen::MatrixXd& P,
                                double h,
                                const std::deque<Eigen::MatrixXd>& prev) {
        return (0.5 * (P + prev[1])) +
               (h / 48.0) * (119 * solver.riccati_equation(P) -
                             99 * solver.riccati_equation(prev[1]) +
                             69 * solver.riccati_equation(prev[2]) -
                             17 * solver.riccati_equation(prev[3]));
    }
};

struct Inglend {
    static Eigen::MatrixXd step(RiccatiSolver& solver, const Eigen::MatrixXd& P,
                                double h, const std::deque<Eigen::MatrixXd>&) {
        Eigen::MatrixXd k1, k2, k3, k4;
        k1.noalias() = solver.riccati_equation(P);
        k2.noalias() = solver.riccati_equation(P + (k1 * (h / 2.0)));
        k3.noalias() = solver.riccati_equation(P + (h / 4.0) * (k1 + k2));
        k4.noalias() = solver.riccati_equation(P + h * (2 * k3 - k2));

        return P + (h / 6.0) * (k1 + 4 * k3 + k4);
    }
};

struct Felberg {
    static Eigen::MatrixXd step(RiccatiSolver& solver, const Eigen::MatrixXd& P,
                                double h, const std::deque<Eigen::MatrixXd>&) {
        Eigen::MatrixXd k1, k2, k3, k4, k5;
        k1.noalias() = solver.riccati_equation(P);
        k2.noalias() = solver.riccati_equation(P + (k1 * (h / 4.0)));
        k3.noalias() = solver.riccati_equation(P + (3 * h / 32.0) * k1 +
                                               (9 * h / 32.0) * k2);
        k4.noalias() = solver.riccati_equation(P + (1932 * h / 2197.0) * k1 -
                                               (7200 * h / 2197.0) * k2 +
                                               (7296 * h / 2197.0) * k3);
        k5.noalias() = solver.riccati_equation(
            P + (439 * h / 216.0) * k1 - (8 * h) * k2 +
            (3680 * h / 513.0) * k3 - (845 * h / 4104.0) * k4);

        return P + h * ((25 * k1 / 216.0) + (1408 * k3 / 2565.0) +
                        (2197 * k4 / 4101.0) - (k5 / 5.0));
    }
};

#endif
