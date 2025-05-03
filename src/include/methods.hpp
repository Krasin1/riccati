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

struct RungeKutta3 {
    static Eigen::MatrixXd step(RiccatiSolver& solver, const Eigen::MatrixXd& P,
                                double h, const std::deque<Eigen::MatrixXd>&) {
        Eigen::MatrixXd k1, k2, k3;
        k1.noalias() = solver.riccati_equation(P);
        k2.noalias() = solver.riccati_equation(P + (0.5 * h) * k1);
        k3.noalias() = solver.riccati_equation(P - h * k1 + 2 * h * k2);
        return P + (h / 6.0) * (k1 + 4 * k2 + k3);
    }
};

struct RungeKutta4 {
    static Eigen::MatrixXd step(RiccatiSolver& solver, const Eigen::MatrixXd& P,
                                double h, const std::deque<Eigen::MatrixXd>&) {
        Eigen::MatrixXd k1, k2, k3, k4;
        k1.noalias() = solver.riccati_equation(P);
        k2.noalias() = solver.riccati_equation(P + (0.5 * h) * k1);
        k3.noalias() = solver.riccati_equation(P + (0.5 * h) * k2);
        k4.noalias() = solver.riccati_equation(P + h * k3);
        return P + (h / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4);
    }
};

struct Adams3 {
    static Eigen::MatrixXd step(RiccatiSolver& solver, const Eigen::MatrixXd& P,
                                double h,
                                const std::deque<Eigen::MatrixXd>& prev) {
        return P + (h / 12.0) * ((23 * solver.riccati_equation(P)) -
                                 (16 * solver.riccati_equation(prev[1])) +
                                 (5 * solver.riccati_equation(prev[2])));
    }
};

struct Adams4 {
    static Eigen::MatrixXd step(RiccatiSolver& solver, const Eigen::MatrixXd& P,
                                double h,
                                const std::deque<Eigen::MatrixXd>& prev) {
        return P + (h / 24.0) * ((55 * solver.riccati_equation(P)) -
                                 (59 * solver.riccati_equation(prev[1])) +
                                 (37 * solver.riccati_equation(prev[2])) -
                                 (9 * solver.riccati_equation(prev[3])));
    }
};

struct Adams5 {
    static Eigen::MatrixXd step(RiccatiSolver& solver, const Eigen::MatrixXd& P,
                                double h,
                                const std::deque<Eigen::MatrixXd>& prev) {
        return P + (h / 720.0) * ((1901 * solver.riccati_equation(P)) -
                                  (2774 * solver.riccati_equation(prev[1])) +
                                  (2616 * solver.riccati_equation(prev[2])) -
                                  (1274 * solver.riccati_equation(prev[3])) +
                                  (251 * solver.riccati_equation(prev[4])));
    }
};

struct Felberg4 {
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

struct Felberg5 {
    static Eigen::MatrixXd step(RiccatiSolver& solver, const Eigen::MatrixXd& P,
                                double h, const std::deque<Eigen::MatrixXd>&) {
        Eigen::MatrixXd k1, k2, k3, k4, k5, k6;
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
        k6.noalias() = solver.riccati_equation(
            P - (8 * h / 27.0) * k1 + 2 * h * k2 - (3544 * h / 2565.0) * k3 +
            (1859 * h / 4101.0) * k4 - (11 * h / 40.0) * k5);
        return P +
               h * ((16 / 135.0) * k1 + (6656 / 12825.0) * k3 +
                    (28561 / 56430.0) * k4 - (9 / 50.0) * k5 + (2 / 55.0) * k6);
    }
};

struct Inglend4 {
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

struct Inglend5 {
    static Eigen::MatrixXd step(RiccatiSolver& solver, const Eigen::MatrixXd& P,
                                double h, const std::deque<Eigen::MatrixXd>&) {
        Eigen::MatrixXd k1, k2, k3, k4, k5, k6;
        k1.noalias() = solver.riccati_equation(P);
        k2.noalias() = solver.riccati_equation(P + (k1 * (h / 2.0)));
        k3.noalias() = solver.riccati_equation(P + (h / 4.0) * (k1 + k2));
        k4.noalias() = solver.riccati_equation(P + h * (2 * k3 - k2));
        k5.noalias() =
            solver.riccati_equation(P + (h / 27.0) * (7 * k1 + 10 * k2 + k4));
        k6.noalias() = solver.riccati_equation(
            P +
            (h / 625.0) * (28 * k1 - 125 * k2 + 546 * k3 + 54 * k4 - 378 * k5));

        return P + (h / 336.0) * (14 * k1 + 35 * k4 + 162 * k5 + 125 * k6);
    }
};

struct Nystrom2 {
    static Eigen::MatrixXd step(RiccatiSolver& solver, const Eigen::MatrixXd& P,
                                double h,
                                const std::deque<Eigen::MatrixXd>& prev) {
        return prev[1] + 2 * h * solver.riccati_equation(P);
    }
};

struct Nystrom3 {
    static Eigen::MatrixXd step(RiccatiSolver& solver, const Eigen::MatrixXd& P,
                                double h,
                                const std::deque<Eigen::MatrixXd>& prev) {
        return prev[1] + (h / 3.0) * (7 * solver.riccati_equation(P) -
                                      2 * solver.riccati_equation(prev[1]) +
                                      solver.riccati_equation(prev[2]));
    }
};

struct Nystrom4 {
    static Eigen::MatrixXd step(RiccatiSolver& solver, const Eigen::MatrixXd& P,
                                double h,
                                const std::deque<Eigen::MatrixXd>& prev) {
        return prev[1] + (h / 3.0) * (8 * solver.riccati_equation(P) -
                                      5 * solver.riccati_equation(prev[1]) +
                                      4 * solver.riccati_equation(prev[2]) -
                                      solver.riccati_equation(prev[3]));
    }
};

struct Milna4 {
    static Eigen::MatrixXd step(RiccatiSolver& solver, const Eigen::MatrixXd& P,
                                double h,
                                const std::deque<Eigen::MatrixXd>& prev) {
        return prev[3] +
               ((4 * h / 3.0) * (2 * solver.riccati_equation(P) -
                                 solver.riccati_equation(prev[1]) +
                                 2 * solver.riccati_equation(prev[2])));
    }
};

struct Milna6 {
    static Eigen::MatrixXd step(RiccatiSolver& solver, const Eigen::MatrixXd& P,
                                double h,
                                const std::deque<Eigen::MatrixXd>& prev) {
        return prev[5] +
               (3 * h / 10.0) * (11 * solver.riccati_equation(P) -
                                 14 * solver.riccati_equation(prev[1]) +
                                 26 * solver.riccati_equation(prev[2]) -
                                 14 * solver.riccati_equation(prev[3]) +
                                 11 * solver.riccati_equation(prev[4]));
    }
};

struct Hemming4 {
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

#endif
