#include <iostream>
#include <math.h>

int main()
{
    double eps = 1;

    while (1.0 + eps > 1.0)
        eps /= 2.0;

    std::cout << "a) eps= " << eps << "\n";

    double x_min = 1;

    while (x_min != x_min + 1)
    {
        x_min *= 2;
    }

    std::cout << "b) x_min= " << x_min << "\n";

    double y_min = 1, num = std::pow(10, 20);

    while (y_min + num != y_min)
    {
        y_min *= 2;
    }

    std::cout << "c) y_min= " << y_min << "\n";

    return 0;
}