#include "include/utils.hpp"

#include <CLI/CLI.hpp>
#include <Eigen/Dense>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

namespace fs = std::filesystem;

// Ввод одной матрицы в файл
Eigen::MatrixXd read_matrix_from_file(const std::string& filepath) {
    std::ifstream file(filepath);
    if (!file.is_open()) {
        throw std::runtime_error("Не удалось открыть файл " + filepath);
    }

    std::vector<std::vector<double>> data;
    std::string line;
    size_t cols = 0;  // Количество столбцов (определяется по первой строке)

    while (std::getline(file, line)) {
        // Пропускаем пустые строки
        if (line.empty()) continue;

        std::istringstream iss(line);
        std::vector<double> row;
        std::string token;

        while (iss >> token) {
            // Проверяем, что число состоит только из цифр, знака и точек 
            bool is_valid = !token.empty();
            bool has_digit = false;
            size_t sign_pos = (token[0] == '-' || token[0] == '+') ? 1 : 0;
            size_t dot_count = 0;

            for (size_t i = sign_pos; i < token.size(); ++i) {
                if (token[i] == '.') {
                    dot_count++;
                    // Не больше одной точки в числе
                    if (dot_count > 1) {
                        is_valid = false;
                        break;
                    }
                } else if (!isdigit(token[i])) {
                    is_valid = false;
                    break;
                } else {
                    has_digit = true;
                }
            }

            if (!is_valid || !has_digit) {
                throw std::runtime_error("Обнаружен лишний символ: " + token);
            }

            try {
                double value = std::stod(token);
                row.push_back(value);
            } catch (...) {
                throw std::runtime_error("Ошибка обработки: " + token);
            }
        }

        if (row.empty()) {
            continue;  // Пропускаем строки без чисел 
        }

        // Проверяем согласованность количества столбцов
        if (cols == 0) {
            cols = row.size();
        } else if (row.size() != cols) {
            throw std::runtime_error("Несовпадение количества столбцов в строке: " + line);
        }

        data.push_back(row);
    }

    file.close();

    if (data.empty()) {
        throw std::runtime_error("Файл " + filepath + " не содержит числовых данных.");
    }

    // Преобразуем в Eigen::MatrixXd
    Eigen::MatrixXd mat(data.size(), cols);
    for (size_t i = 0; i < data.size(); ++i) {
        for (size_t j = 0; j < data[i].size(); ++j) {
            mat(static_cast<long>(i), static_cast<long>(j)) = data[i][j];
        }
    }

    return mat;
}

// Вывод прогресса раз в 2 секунды
void progress(Config cfg, double error, int step, double t,
              Eigen::MatrixXd& P) {
    static auto start_time = std::chrono::system_clock::now();
    auto current_time = std::chrono::system_clock::now();
    auto elapsed_time = std::chrono::duration_cast<std::chrono::seconds>(
        current_time - start_time);

    if (cfg.manual || elapsed_time.count() >= 2) {
        std::cout << "\033[H\033[2J";
        std::cout << "Метод: " << cfg.method
                  << " | Требуемая ошибка: " << cfg.target_error
                  << " | Ошибка: " << error << "\t\r\nh: " << cfg.h
                  << " | Шаг : " << step << " | t : " << t << " из "
                  << cfg.t_max << "\t\t\n";

        std::cout << "P(первые 8x8 элементов):\n";
        int rows_to_show = std::min(static_cast<int>(P.rows()), 8);
        int cols_to_show = std::min(static_cast<int>(P.cols()), 8);

        for (int i = 0; i < rows_to_show; ++i) {
            for (int j = 0; j < cols_to_show; ++j) {
                std::cout << std::setw(10) << std::setprecision(4) << P(i, j)
                          << " ";
            }
            std::cout << "\n";
        }
        start_time = current_time;
    }

    int ch;
    while (cfg.manual) {
        ch = std::cin.get();
        if (ch == '\n') break;
    }
}

// записывает результаты в файл и выводит в терминал
void show_results(Config cfg, double error, double steps,
                  std::chrono::milliseconds elapsed_time, Eigen::MatrixXd& P) {
    fs::path folder_path = "results";
    fs::path result_file_path = folder_path / "output.txt";

    if (!fs::exists(folder_path)) fs::create_directory(folder_path);
    std::ofstream out(result_file_path);

    if (out.is_open()) {
        out << cfg.method << " ; t = " << cfg.t0 << " ; t_end = " << cfg.t_max
            << " ; h = " << cfg.h << " ; error_value = " << cfg.target_error
            << "\n"
            << "error = " << error << " ; last_step = " << steps
            << "\ntime: " << elapsed_time.count() / 1000 << " sec. "
            << elapsed_time.count() % 1000 << " ms."
            << "\n\n";
        out << std::fixed << std::setprecision(7);
        out << P << "\n";
        out.close();
    } else {
        std::cout << "Unable to open file\n";
    }

    std::cout << "\033[2J\033[H\nИсходные данные: ";
    std::cout << cfg.method << "; t = " << cfg.t0 << " ; t_end = " << cfg.t_max
              << " ; h = " << cfg.h << " ; error_value = " << cfg.target_error
              << "\n";
    std::cout << "Последняя ошибка = " << error << " ; Шаг = " << steps << "\n";
    std::cout << "\033[32mMатрица P записана в output.txt\033[0m\n";
    std::cout << "Время работы программы: " << elapsed_time.count() / 1000
              << " сек. " << elapsed_time.count() % 1000 << " мс.\n";
}

// Рисует график если флаг draw == true
void draw_graph(std::vector<double>* error, Config cfg) {
    FILE* gnuplot = popen("gnuplot -persist", "w");
    if (!gnuplot) {
        std::cerr << "gnuplot не удалось запустить\n";
        return;
    }

    fs::path folder = "results/png";
    if (!fs::exists(folder)) fs::create_directory(folder);

    fprintf(gnuplot, "set terminal pngcairo\n");
    fprintf(gnuplot, "set output 'results/png/%s.png'\n", cfg.method.c_str());
    fprintf(gnuplot, "set title \"Ошибка на каждом шаге\"\n");
    fprintf(gnuplot, "set xlabel 'Шаги'\n");
    fprintf(gnuplot, "set ylabel 'Ошибка'\n");
    fprintf(gnuplot, "plot '-' with lines title 'error'\n");

    for (size_t i = 0; i < error->size(); i++) {
        fprintf(gnuplot, "%ld %f\n", i, error->at(i));
    }
    fprintf(gnuplot, "e\n");
    fprintf(gnuplot, "set output\n");
    pclose(gnuplot);
}

void check_nan(const Eigen::MatrixXd& P, int step, double t) {
    if (std::isnan(P(0, 0))) {
        std::string message =
            "Матрица P заполнилась -nan\nШаг : " + std::to_string(step) +
            " | t : " + std::to_string(t) + "\n";
        throw std::runtime_error(message);
    }
}

Config parse_cli(int argc, char** argv) {
    Config opts;
    CLI::App app{"Riccati Equation Solver"};

    app.add_option("--method", opts.method,
                   "Метод решения (runge3, runge4,"
                   "adams3, adams4, adams5, "
                   "milna4, milna6, "
                   "nystrom2, nystrom3, nystrom4, "
                   "felberg4, felberg5, "
                   "inglend4, inglend5, "
                   "hemming4)")
        ->required();
    app.add_option("--h", opts.h, "Шаг интегрирования >= 0");
    app.add_option("--t0", opts.t0, "Начальное время (по умолчанию 0)");
    app.add_option("--t_max", opts.t_max, "Конечное время (по умолчанию 10)");
    app.add_option("--error", opts.target_error,
                   "Требуемая погрешность (по умолчанию 0.001)");
    app.add_option("--max_steps", opts.max_steps,
                   "Максимум шагов (по умолчанию 200)");
    app.add_option("--threads", opts.threads,
                   "Количество потоков (по умолчанию 1)\n"
                   "Для небольших матриц большое кол-во потоков негативно "
                   "повлияют на скорость.");
    app.add_flag("--draw", opts.draw,
                 "Рисовать график ошибки (использует gnuplot)");
    app.add_flag("--manual", opts.manual, "Пошаговое выполнение вычислений");
    app.add_flag("--step_time", opts.step_time,
                 "Отображать время выполнения одного шага");

    if (argc < 2) {
        std::cout << app.help() << "\n";
        exit(0);
    }

    try {
        app.parse(argc, argv);
    } catch (const CLI::ParseError& e) {
        app.exit(e);
        exit(1);
    }

    return opts;
}
