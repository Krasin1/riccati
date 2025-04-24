#include <iostream>

#include "include/fabric.hpp"
#include "include/utils.hpp"

int main(int argc, char* argv[]) {
    try {
        // Параметры интегрирования из параметров
        Config cfg = parse_cli(argc, argv);

        // включаем многопоток
        Eigen::setNbThreads(cfg.threads);

        Eigen::MatrixXd E = read_matrix_from_file("data/E.dat");
        Eigen::MatrixXd A = read_matrix_from_file("data/A.dat");
        Eigen::MatrixXd B = read_matrix_from_file("data/B.dat");
        Eigen::MatrixXd Q = read_matrix_from_file("data/Q.dat");
        Eigen::MatrixXd initial_P = Eigen::MatrixXd::Zero(E.rows(), E.cols());

        // Список методов
        //  - runge
        //  - adams
        //  - milna
        //  - nystrom
        //  - hemming
        //  - inglend
        //  - felberg
        auto solver = create_solver(cfg.method, E, A, B, Q, initial_P);

        // Замеры времени
        auto begin = std::chrono::system_clock::now();

        // Решаем
        Result result = solver->solve(cfg);

        auto end = std::chrono::system_clock::now();
        auto duration =
            std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);

        // Подставляем найденную матрицу в уравнение риккати -> ответ записываем
        // в файл
        solver->verify_solution(result.P);

        // результаты находятся в папке results
        show_results(cfg, result.last_error, result.step, duration, result.P);

        // Рисуем график если передали агрумент draw при запуске программы
        if (cfg.draw) draw_graph(result.points);
    } catch (std::exception& e) {
        std::cout << e.what() << '\n';
    }

    return 0;
}
