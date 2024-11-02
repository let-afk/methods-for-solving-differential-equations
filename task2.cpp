#include <iostream>
#include <vector>
#include <math.h>

double recurrent_integral(double prevValue, int n)
{
    const double nextValue = 1.0 / double(n) - 7.0 * prevValue;
    return nextValue;
}

double reverse_recurrent_integral(double prevValue, int n)
{
    const double nextValue = 1.0 / (7.0 * n) - prevValue / 7.0;
    return nextValue;
}

double specifiedFunction(double x, int n)
{
    return std::pow(x, n) / (x + 7);
}

double partialAmountSpecificAmount(int countSteps, int n)
{
    double partSum = 0;
    for (int step = 0; step < countSteps; ++step)
    {
        partSum += (1.0 / countSteps) * specifiedFunction((double(step) / countSteps) + 1.0 / (countSteps * 2.0), n);
    }

    return partSum;
}

int main()
{
    double initial_value = std::log(8.0 / 7.0);
    double recurrent_value = initial_value;

    for (int n = 1; n <= 31; ++n)
    {
        recurrent_value = recurrent_integral(recurrent_value, n);
    }

    std::cout << "a) I31= " << recurrent_value << "\n";

    double reverse_recurrent_value = 0;

    for (int n = 61; n >= 32; --n)
    {
        reverse_recurrent_value = reverse_recurrent_integral(reverse_recurrent_value, n);
    }

    std::cout << "b) I31= " << reverse_recurrent_value << "\n";

    double sum_value = partialAmountSpecificAmount(1000, 31);

    std::cout << "c) I31= " << sum_value << "\n";

    return 0;
}