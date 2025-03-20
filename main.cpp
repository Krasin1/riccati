#include <iostream>
#include <thread>

#include "functions/func.h"
#include "methods/methods.h"

bool draw = false;
bool manual = false;

int main(int argc, char* argv[]) {
    try {
        set_flags(argc, argv);

        // включаем многопоток
        Eigen::setNbThreads(std::thread::hardware_concurrency());

        Eigen::MatrixXd E = read_matrix_from_file("!mat/E.dat");
        Eigen::MatrixXd A = read_matrix_from_file("!mat/A.dat");
        Eigen::MatrixXd B = read_matrix_from_file("!mat/B.dat");
        Eigen::MatrixXd Q = read_matrix_from_file("!mat/Q.dat");
        Eigen::MatrixXd initial_P = Eigen::MatrixXd::Zero(E.rows(), E.cols());

        // Параметры интегрирования
        double t0 = 0.0;
        double t_max = 10.0;
        double h = 0.1;
        double target_error = 0.001;

        // ввод из терминала
        input_data(t0, t_max, h, target_error);

        // Создаём решатель
        FelbergSolver solver(E, A, B, Q, initial_P);

        // Замеры времени
        auto begin = std::chrono::system_clock::now();

        // Решаем
        SolverResult result = solver.solve(t0, t_max, h, target_error);

        auto end = std::chrono::system_clock::now();
        auto duration =
            std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);

        // Подставляем найденную матрицу в уравнение риккати -> ответ записываем
        // в файл
        solver.verify_solution(result.P);

        // результаты находятся в папке Results
        show_results(t0, t_max, h, target_error, result.last_error, result.step,
                     duration, result.P);

        // Рисуем график если передали агрумент draw при запуске программы
        if (draw) draw_graph(result.points);

    } catch (std::exception& e) {
        std::cout << e.what() << '\n';
    }

    return 0;
}
