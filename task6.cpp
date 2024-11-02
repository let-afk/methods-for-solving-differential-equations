#include <iostream>
#include <fstream>
#include <numbers>
#include <math.h>
#include <vector>

std::vector<double> specified_function(std::vector<double> y_vector)
{
    std::vector<double> temp;
    temp.push_back(y_vector[1]);
    temp.push_back(-y_vector[0]);
    return temp;
}

void runge_kutta_method(double h, std::vector<double> y_vector, std::string file_name)
{
    std::ofstream out;
    out.open(file_name);
    double t = 0, t_max = 100 * M_PI;
    const int s = 7;

    std::vector<std::vector<double>> butcher_table_a = {{0.2, 0, 0, 0, 0, 0},
                                                        {3.0 / 40, 9.0 / 40, 0, 0, 0, 0},
                                                        {44.0 / 45, -56.0 / 15, 32.0 / 9, 0, 0, 0},
                                                        {19372.0 / 6561, -25360.0 / 2187, 64448.0 / 6561, -212.0 / 729, 0, 0},
                                                        {9017.0 / 3168, -355.0 / 33, 46732.0 / 5247, 49.0 / 176, -5103.0 / 18656, 0},
                                                        {35.0 / 384, 0, 500.0 / 1113, 125.0 / 192, -2187.0 / 6784, 11.0 / 84}};
    std::vector<double> butcher_table_b = {5179.0 / 57600, 0, 7571.0 / 16695, 393.0 / 640, -92097.0 / 339200, 187.0 / 2100, 1.0 / 40};
    std::vector<double> butcher_table_c = {0, 0.2, 0.3, 0.8, 8.0 / 9, 1, 1};
    std::vector<std::vector<double>> butcher_table_k(7);

    out << "x;z\n";

    while (t < t_max)
    {
        out << y_vector[0] << ";" << y_vector[1] << "\n";

        butcher_table_k[0] = specified_function(y_vector);

        for (int i = 0; i < s; ++i)
        {
            std::vector<double> sum = y_vector;

            for (int j = 0; j < i; ++j)
            {
                sum[0] += h * butcher_table_a[i - 1][j] * butcher_table_k[j][0];
                sum[1] += h * butcher_table_a[i - 1][j] * butcher_table_k[j][1];
            }

            butcher_table_k[i] = specified_function(sum);
        }

        std::vector<double> y_additional_sum(2);
        for (int i = 0; i < s; i++)
        {
            y_additional_sum[0] += butcher_table_b[i] * butcher_table_k[i][0];
            y_additional_sum[1] += butcher_table_b[i] * butcher_table_k[i][1];
        }

        y_vector[0] += h * y_additional_sum[0];
        y_vector[1] += h * y_additional_sum[1];

        t += h;
    }

    out.close();
}

int main()
{
    std::vector<double> y_vector = {0, 1.0};
    runge_kutta_method(1.0, y_vector, "task6-h1.csv");
    runge_kutta_method(0.1, y_vector, "task6-h01.csv");
    runge_kutta_method(0.01, y_vector, "task6-h001.csv");

    return 0;
}