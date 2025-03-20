#ifndef FUNC_H
#define FUNC_H

#include <chrono>
#include <Eigen/Dense>
#include <string>

// Функция для чтения матрицы из файла
Eigen::MatrixXd read_matrix_from_file(const std::string& filepath);

// Функция для отображения прогресса в терминал
void progress(double target_error, double error, double t, double h, int step,
              double t_max, Eigen::MatrixXd& P);

// Ввод числа с обработкой
template <class T>
T input_number(const std::string& message, const std::string& letter);

// ввод параметров интегрирования из терминала
void input_data(double& t_0, double& t_max, double& h, double& error);

// Параметры программы при запуске
void set_flags(int argc, char* argv[]);

// записывает результаты в файл и выводит в терминал
void show_results(double t0, double t_max, double h, double target_error,
                  double error, double steps,
                  std::chrono::milliseconds elapsed_time, Eigen::MatrixXd& P);

// Рисует график если флаг draw == true
void draw_graph(std::vector<double>* error);

#endif  // !FUNC_H
