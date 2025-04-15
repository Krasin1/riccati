#ifndef SOLVER_H
#define SOLVER_H

#include <Eigen/Dense>
#include <deque>

// Структура для хранения результата
struct SolverResult {
    Eigen::MatrixXd P;
    int step;
    double last_error;
    std::vector<double>* points;
};

// Абстрактный базовый класс для решателей уравнения Риккати
class RiccatiSolver {
   public:
    Eigen::MatrixXd E, A, BRB, Q, initial_P;

    RiccatiSolver(Eigen::MatrixXd& E, Eigen::MatrixXd& A, Eigen::MatrixXd& B,
                  Eigen::MatrixXd& Q, Eigen::MatrixXd& initial_P) {
        this->E = E;
        this->A = A;
        this->BRB = B * B.transpose();
        this->Q = Q;
        this->initial_P = initial_P;
    }

    // виртуальная функция, которую все методы реализуют для решения уравнения
    virtual SolverResult solve(double t0, double t_max, double h,
                               double target_error, int max_steps) = 0;
    // Наше матричное уравнение Риккати
    virtual Eigen::MatrixXd riccati_equation(const Eigen::MatrixXd& P);
    // Найденную матрицу подставляем в уравнение Риккати и выводим в файл
    virtual Eigen::MatrixXd verify_solution(const Eigen::MatrixXd& P);
    // С помощью метода Рунге-Кутты получаем начальные точки для некоторых
    // методов
    virtual std::deque<Eigen::MatrixXd> acceleration_points(int count, double h, Eigen::MatrixXd& P_initial);
    virtual std::deque<Eigen::MatrixXd> acceleration_points(int count, double h);

    virtual ~RiccatiSolver() {}
};

#endif  // SOLVER_H
