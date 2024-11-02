#include <iostream>
#include <fstream>
#include <numbers>
#include <math.h>

void get_info_for_draw(std::string file_name, double h, double x0, double z0)
{
    std::ofstream out;
    out.open(file_name);
    double current_x = x0, current_z = z0, t = 0.;

    out << "x;z\n";

    while (t < 100 * M_PI)
    {
        out << current_x << ";" << current_z << "\n";
        const double prev_x = current_x;
        current_x += h * current_z;
        current_z -= h * prev_x;
        t += h;
    }

    out.close();
}

int main()
{
    get_info_for_draw("task5-h1.csv", 1., 0.0, 1.0);
    get_info_for_draw("task5-h01.csv", 0.1, 0.0, 1.0);
    get_info_for_draw("task5-h001.csv", 0.01, 0.0, 1.0);

    return 0;
}