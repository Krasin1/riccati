#ifndef UTILS_HPP
#define UTILS_HPP

#include <Eigen/Dense>
#include <chrono>
#include <string>
#include <vector>

struct Config {
    double t0 = 0.0;
    double t_max = 10.0;
    double h = 0.01;
    double target_error = 0.001;
    int max_steps = 200;
    int threads = 1;
    bool draw = false;
    bool manual = false;
    bool step_time = false;
    std::string method = "runge4";
};

// Структура для хранения результата
struct Result {
    Eigen::MatrixXd P;
    int step;
    double last_error;
    std::vector<double>* points;
};

Config parse_cli(int argc, char* argv[]);

// Чтение матрицы из файла
Eigen::MatrixXd read_matrix_from_file(const std::string& filepath);

// Вывод прогресса вычислений
void progress(Config cfg, double error, int step, double t, Eigen::MatrixXd& P);

// Вывод результатов вычислений
void show_results(Config cfg, double error, double steps,
                  std::chrono::milliseconds elapsed_time, Eigen::MatrixXd& P);

// Построение графика ошибки (если разрешено)
void draw_graph(std::vector<double>* error, Config cfg);

// Проверка на NaN в матрице
void check_nan(const Eigen::MatrixXd& P, int step, double t);

#endif
