#include <iostream>
#include <fstream>
#include <math.h>

double sub_func_a(double t)
{
    return 2 * t;
}

double sub_func_b(double y)
{
    return y;
}

double sub_func_c(double t)
{
    return 3 * t * t + std::cos(t);
}

void get_info_for_draw(std::string file_name, double h, double (*func)(double), double y0)
{
    std::ofstream out;
    out.open(file_name);
    const int END = 10;
    double current_y = y0, t = 0.;

    out << "t;y\n";

    for (int step = 0; step < END * (1 / h); ++step)
    {
        out << t << ";" << current_y << "\n";
        current_y += h * func(t);
        t += h;
    }

    std::cout << std::abs(current_y - (t * t));

    out.close();
}

void get_info_for_draw_dependent_on_y(std::string file_name, double h, double (*func)(double), double y0)
{
    std::ofstream out;
    out.open(file_name);
    const int END = 10;
    double current_y = y0, t = 0.;

    out << "t;y\n";

    for (int step = 0; step < END * (1 / h); ++step)
    {
        out << t << ";" << current_y << "\n";
        current_y += h * func(current_y);
        t += h;
    }

    out.close();
}

int main()
{
    get_info_for_draw("task4-h1.csv", 1., sub_func_a, 0.0);
    get_info_for_draw("task4-h01.csv", 0.1, sub_func_a, 0.0);
    get_info_for_draw("task4-h001.csv", 0.01, sub_func_a, 0.0);

    get_info_for_draw_dependent_on_y("task4-b-h1.csv", 1., sub_func_b, 1.0);
    get_info_for_draw_dependent_on_y("task4-b-h01.csv", 0.1, sub_func_b, 1.0);
    get_info_for_draw_dependent_on_y("task4-b-h001.csv", 0.01, sub_func_b, 1.0);

    get_info_for_draw("task4-c-h1.csv", 1., sub_func_c, 5.0);
    get_info_for_draw("task4-c-h01.csv", 0.1, sub_func_c, 5.0);
    get_info_for_draw("task4-c-h001.csv", 0.01, sub_func_c, 5.0);

    return 0;
}