#ifndef FABRIC_HPP
#define FABRIC_HPP

#include <Eigen/Dense>
#include <memory>

#include "methods.hpp"

inline std::unique_ptr<RiccatiSolver> create_solver(
    std::string& name, Eigen::MatrixXd& E, Eigen::MatrixXd& A,
    Eigen::MatrixXd& B, Eigen::MatrixXd& Q, Eigen::MatrixXd& P0) {
    if (name == "runge3")
        return std::make_unique<Solver<RungeKutta3>>(E, A, B, Q, P0);
    if (name == "runge4")
        return std::make_unique<Solver<RungeKutta4>>(E, A, B, Q, P0);

    if (name == "adams3")
        return std::make_unique<Solver<Adams3>>(E, A, B, Q, P0);
    if (name == "adams4")
        return std::make_unique<Solver<Adams4>>(E, A, B, Q, P0);
    if (name == "adams5")
        return std::make_unique<Solver<Adams5>>(E, A, B, Q, P0);

    if (name == "milna4")
        return std::make_unique<Solver<Milna4>>(E, A, B, Q, P0);
    if (name == "milna6")
        return std::make_unique<Solver<Milna6>>(E, A, B, Q, P0);

    if (name == "nystrom2")
        return std::make_unique<Solver<Nystrom2>>(E, A, B, Q, P0);
    if (name == "nystrom3")
        return std::make_unique<Solver<Nystrom3>>(E, A, B, Q, P0);
    if (name == "nystrom4")
        return std::make_unique<Solver<Nystrom4>>(E, A, B, Q, P0);

    if (name == "hemming4")
        return std::make_unique<Solver<Hemming4>>(E, A, B, Q, P0);

    if (name == "felberg4")
        return std::make_unique<Solver<Felberg4>>(E, A, B, Q, P0);
    if (name == "felberg5")
        return std::make_unique<Solver<Felberg5>>(E, A, B, Q, P0);

    if (name == "inglend4")
        return std::make_unique<Solver<Inglend4>>(E, A, B, Q, P0);
    if (name == "inglend5")
        return std::make_unique<Solver<Inglend5>>(E, A, B, Q, P0);

    throw std::runtime_error("Неизвестный метод интегрирования: " + name);
}

#endif
