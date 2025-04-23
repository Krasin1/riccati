#include "func.hpp"

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>

namespace fs = std::filesystem;

extern bool draw;
extern bool manual;

// Ввод одной матрицы в файл
Eigen::MatrixXd read_matrix_from_file(const std::string& filepath) {
    std::ifstream file(filepath);
    if (!file.is_open()) {
        throw std::runtime_error("Не удалось открыть файл " + filepath);
    }
    std::vector<std::vector<double>> data;
    std::string line;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::vector<double> row;
        double value;
        while (iss >> value) row.push_back(value);
        if (!row.empty()) data.push_back(row);
    }
    file.close();

    if (data.empty()) {
        throw std::runtime_error("Файл " + filepath +
                                 " пустой или содержит неверные данные.");
    }

    Eigen::MatrixXd mat(data.size(), data[0].size());
    for (size_t i = 0; i < data.size(); ++i) {
        for (size_t j = 0; j < data[i].size(); ++j) {
            mat((long)i, (long)j) = data[i][j];
        }
    }
    return mat;
}
// Вывод прогресса раз в 2 секунды
void progress(double target_error, double error, double t, double h, int step,
              double t_max, Eigen::MatrixXd& P) {
    static auto start_time = std::chrono::system_clock::now();
    auto current_time = std::chrono::system_clock::now();
    auto elapsed_time = std::chrono::duration_cast<std::chrono::seconds>(
        current_time - start_time);

    if (manual || elapsed_time.count() >= 2) {
        std::cout << "\033[H\033[2J";
        std::cout << "Требуемая ошибка: " << target_error
                  << " | Ошибка: " << error << "\t\r\nh: " << h
                  << " | Шаг : " << step << " | t : " << t << " из " << t_max
                  << "\t\t\n";

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
    while (manual) {
        ch = std::cin.get();
        if (ch == '\n') break;
    }
}

// Ввод числа с обработкой
template <class T>
T input_number(const std::string& message, const std::string& letter) {
    T a;
    while (true) {
        std::cout << message << " " << letter << ": ";
        if (std::cin >> a && a >= 0) {
            break;
        } else {
            std::cin.clear();
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        }
    }
    return a;
}

// ввод параметров интегрирования из терминала
void input_data(double& t_0, double& t_max, double& h, double& error) {
    // t_0 = input_number<double>("\033[2J\033[HВведите начальное время", "t");
    // t_max =
    //     input_number<double>("\033[2J\033[HВведите конечное время", "t_max");
    h = input_number<double>("\033[2J\033[HВведите шаг", "h");
    if (h <= 0) {
        std::cout << "Шаг не может быть <= 0\n";
        exit(1);
    }
    error = input_number<double>("\033[2J\033[HВведите необходимую ошибку", "");
    if (error <= 0) {
        std::cout << "Точность должна быть больше нуля\n";
        exit(1);
    }
}

// записывает результаты в файл и выводит в терминал
void show_results(double t_0, double t_max, double h, double target_error,
                  double error, double steps,
                  std::chrono::milliseconds elapsed_time, Eigen::MatrixXd& P) {
    fs::path folder_path = "results";
    fs::path result_file_path = folder_path / "output.txt";

    if (!fs::exists(folder_path)) fs::create_directory(folder_path);
    std::ofstream out(result_file_path);

    if (out.is_open()) {
        out << "t = " << t_0 << " ; t_end = " << t_max << " ; h = " << h
            << " ; error_value = " << target_error << "\n"
            << "error = " << error << " ; last_step = " << steps
            << "\ntime: " << elapsed_time.count() / 1000 << " sec. "
            << elapsed_time.count() % 1000 << " ms."
            << "\n\n";
        out << P << "\n";
        out.close();
    } else {
        std::cout << "Unable to open file\n";
    }

    std::cout << "\033[2J\033[H\nИсходные данные: ";
    std::cout << "t = " << t_0 << " ; t_end = " << t_max << " ; h = " << h
              << " ; error_value = " << target_error << "\n";
    std::cout << "Последняя ошибка = " << error << " ; Шаг = " << steps << "\n";
    std::cout << "\033[32mMатрица P записана в output.txt\033[0m\n";
    std::cout << "Время работы программы: " << elapsed_time.count() / 1000
              << " сек. " << elapsed_time.count() % 1000 << " мс.\n";
}

// Рисует график если флаг draw == true
void draw_graph(std::vector<double>* error) {
    FILE* gnuplot = popen("gnuplot -persist", "w");
    if (!gnuplot) {
        std::cerr << "gnuplot не удалось запустить\n";
        return;
    }
    fprintf(gnuplot, "set title \"Ошибка на каждом шаге\"\n");
    fprintf(gnuplot, "set xlabel 'Шаги'\n");
    fprintf(gnuplot, "set ylabel 'Ошибка'\n");
    fprintf(gnuplot, "plot '-' with lines title 'error'\n");

    for (size_t i = 0; i < error->size(); i++) {
        fprintf(gnuplot, "%ld %f\n", i, error->at(i));
    }
    fprintf(gnuplot, "e\n");
    pclose(gnuplot);
}

// Параметры программы при запуске
void set_flags(int argc, char* argv[]) {
    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "draw") {
            draw = true;
        } else if (arg == "manual") {
            manual = true;
        }
    }
}

void check_nan(const Eigen::MatrixXd& P, int step, double t) {
    if (std::isnan(P(0, 0))) {
        std::string message =
            "Матрица P заполнилась -nan\nШаг : " + std::to_string(step) +
            " | t : " + std::to_string(t) + "\n";
        throw std::runtime_error(message);
    }
}
