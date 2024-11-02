#include <iostream>
#include <stdio.h>
#include <math.h>

double specifiedFunction(double x)
{
    return 5 + x * x * x + std::sin(x);
}

int main()
{
    double x0 = 0.;
    double h = 1.;
    double exact_value_derivative = 1.0;

    std::cout << "h:                        R1: \n";

    for (int step = 0; step <= 20; ++step)
    {
        double r1 = exact_value_derivative - (specifiedFunction(x0 + h) - specifiedFunction(x0)) / h;

        std::printf("%.20f    %.20f \n", h, r1);

        h /= 10.;
    }

    return 0;
}