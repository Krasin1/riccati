#ifndef SOLVER_HPP
#define SOLVER_HPP

#include <Eigen/Dense>
#include <deque>

#include "utils.hpp"

// Абстрактный базовый класс для решателей уравнения Риккати
class RiccatiSolver {
   public:
    RiccatiSolver(Eigen::MatrixXd& E, Eigen::MatrixXd& A, Eigen::MatrixXd& B,
                  Eigen::MatrixXd& Q, Eigen::MatrixXd& initial_P) {
        this->E_ = E;
        this->A_ = A;
        this->BRB_ = B * B.transpose();
        this->Q_ = Q;
        this->initial_P_ = initial_P;
    }
    Result solve(Config cfg);

    int get_acceleration_points(const std::string& method);
    // Наше матричное уравнение Риккати
    virtual Eigen::MatrixXd riccati_equation(const Eigen::MatrixXd& P);
    // Найденную матрицу подставляем в уравнение Риккати и выводим в файл
    virtual Eigen::MatrixXd verify_solution(const Eigen::MatrixXd& P);
    // С помощью метода Рунге-Кутты получаем начальные точки для некоторых
    // методов
    virtual std::deque<Eigen::MatrixXd> acceleration_points(
        int count, double h, Eigen::MatrixXd& P_initial);

    virtual std::deque<Eigen::MatrixXd> acceleration_points(int count,
                                                            double h);

    virtual ~RiccatiSolver() {}

   protected:
    Eigen::MatrixXd E_, A_, BRB_, Q_, initial_P_;

    // виртуальная функция, которую все методы реализуют для решения уравнения
    virtual Eigen::MatrixXd update_step(
        const Eigen::MatrixXd& P, double h,
        const std::deque<Eigen::MatrixXd>& prev_points) = 0;
};

#endif
