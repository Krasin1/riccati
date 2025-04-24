#ifndef FABRIC_HPP
#define FABRIC_HPP

#include <Eigen/Dense>
#include <memory>

#include "methods.hpp"

inline std::unique_ptr<RiccatiSolver> create_solver(
    std::string& name, Eigen::MatrixXd& E, Eigen::MatrixXd& A,
    Eigen::MatrixXd& B, Eigen::MatrixXd& Q, Eigen::MatrixXd& P0) {
    if (name == "runge") {
        return std::make_unique<Solver<RungeKutta>>(E, A, B, Q, P0);
    }
    if (name == "adams") {
        return std::make_unique<Solver<Adams>>(E, A, B, Q, P0);
    }
    if (name == "milna") {
        return std::make_unique<Solver<Milna>>(E, A, B, Q, P0);
    }
    if (name == "nystrom") {
        return std::make_unique<Solver<Nystrom>>(E, A, B, Q, P0);
    }
    if (name == "hemming") {
        return std::make_unique<Solver<Hemming>>(E, A, B, Q, P0);
    }
    if (name == "inglend") {
        return std::make_unique<Solver<Inglend>>(E, A, B, Q, P0);
    }
    if (name == "felberg") {
        return std::make_unique<Solver<Felberg>>(E, A, B, Q, P0);
    }
    throw std::runtime_error("Неизвестный метод интегрирования");
}

#endif
