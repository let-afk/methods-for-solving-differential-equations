#include <iostream>
#include <cmath>
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <fstream>

using namespace std;

double dx_dt(double t, double x, double z)
{
    return z;
}

double dz_dt(double t, double x, double z)
{
    return (1 - (4.0 / 3) * x * (1 + x)) * z - x + 5 * std::cos(2 * t);
}

void runge_kutta(double x, double z, double h, double t, double array_Y[4])
{
    double k1x, k1z, k2x, k2z, k3x, k3z, k4x, k4z, k5x, k5z, k6x, k6z, k7x, k7z;

    k1x = dx_dt(t, x, z);
    k1z = dz_dt(t, x, z);
    k2x = dx_dt(t + (1.0 / 5) * h, x + h * (1.0 / 5) * k1x, z + h * (1.0 / 5) * k1z);
    k2z = dz_dt(t + (1.0 / 5) * h, x + h * (1.0 / 5) * k1x, z + h * (1.0 / 5) * k1z);
    k3x = dx_dt(t + (3. / 10) * h, x + h * (3. / 40) * k1x + h * (9. / 40) * k2x, z + h * (3. / 40) * k1z + h * (9. / 40) * k2z);
    k3z = dz_dt(t + (3. / 10) * h, x + h * (3. / 40) * k1x + h * (9. / 40) * k2x, z + h * (3. / 40) * k1z + h * (9. / 40) * k2z);

    k4x = dx_dt(t + (4. / 5) * h, x + h * (44. / 55) * k1x + h * (-1 * 56. / 15) * k2x + h * (32. / 9) * k3x, z + h * (44. / 55) * k1z + h * (-1 * 56. / 15) * k2z + h * (32. / 9) * k3z);
    k4z = dz_dt(t + (4. / 5) * h, x + h * (44. / 55) * k1x + h * (-1 * 56. / 15) * k2x + h * (32. / 9) * k3x, z + h * (44. / 55) * k1z + h * (-1 * 56. / 15) * k2z + h * (32. / 9) * k3z);

    k5x = dx_dt(t + (8. / 9) * h, x + h * ((19372. / 6561) * k1x + (-1 * 25360. / 2187) * k2x + (64448. / 6561) * k3x + (-1 * 212. / 729) * k4x),
                z + h * ((19372. / 6561) * k1z + (-25360. / 2187) * k2z + (64448. / 6561) * k3z + (-212. / 729) * k4z));
    k5z = dz_dt(t + (8. / 9) * h, x + h * ((19372. / 6561) * k1x + (-25360. / 2187) * k2x + (64448. / 6561) * k3x + (-1 * 212. / 729) * k4x),
                z + h * ((19372. / 6561) * k1z + (-1 * 25360. / 2187) * k2z + (64448. / 6561) * k3z + (-212. / 729) * k4z));

    k6x = dx_dt(t + h, x + h * ((9017. / 3168) * k1x + (-1 * 355. / 33) * k2x + (46732. / 5247) * k3x + (49. / 176) * k4x + (-5103. / 18656) * k5x),
                z + h * ((9017. / 3168) * k1z + (-355. / 33) * k2z + (46732. / 5247) * k3z + (49. / 176) * k4z + (-1 * 5103. / 18656) * k5z));
    k6z = dz_dt(t + h, x + h * ((9017. / 3168) * k1x + (-1 * 355. / 33) * k2x + (46732. / 5247) * k3x + (49. / 176) * k4x + (-5103. / 18656) * k5x),
                z + h * ((9017. / 3168) * k1z + (-355. / 33) * k2z + (46732. / 5247) * k3z + (49. / 176) * k4z + (-1 * 5103. / 18656) * k5z));

    k7x = dx_dt(t + h, x + h * ((35. / 384) * k1x + (0.) * k2x + (500. / 1113) * k3x + (125. / 192) * k4x + (-1 * 2187. / 6784) * k5x + (11. / 84) * k6x),
                z + h * ((35. / 384) * k1z + (0.) * k2z + (500. / 1113) * k3z + (125. / 192) * k4z + (-1 * 2187. / 6784) * k5z + (11. / 84) * k6z));
    k7z = dz_dt(t + h, x + h * ((35. / 384) * k1x + (0.) * k2x + (500. / 1113) * k3x + (125. / 192) * k4x + (-1 * 2187. / 6784) * k5x + (11. / 84) * k6x),
                z + h * ((35. / 384) * k1z + (0.) * k2z + (500. / 1113) * k3z + (125. / 192) * k4z + (-1 * 2187. / 6784) * k5z + (11. / 84) * k6z));

    array_Y[0] = x + h * ((35.0 / 384) * k1x + (0.0) * k2x + (500.0 / 1113) * k3x + (125.0 / 192) * k4x + (-2187.0 / 6784) * k5x + (11.0 / 84) * k6x + (0.0) * k7x);
    array_Y[1] = z + h * ((35.0 / 384) * k1z + (0.0) * k2z + (500.0 / 1113) * k3z + (125.0 / 192) * k4z + (-2187.0 / 6784) * k5z + (11.0 / 84) * k6z + (0.0) * k7z);
    array_Y[2] = x + h * ((5179.0 / 57600) * k1x + (0.0) * k2x + (7571.0 / 16695) * k3x + (393.0 / 640) * k4x + (-92097.0 / 339200) * k5x + (187.0 / 2100) * k6x + (1.0 / 40) * k7x);
    array_Y[3] = z + h * ((5179.0 / 57600) * k1z + (0.0) * k2z + (7571.0 / 16695) * k3z + (393.0 / 640) * k4z + (-92097.0 / 339200) * k5z + (187.0 / 2100) * k6z + (1.0 / 40) * k7z);
}

double def_err(const double *array_Y)
{
    return max(fabs(array_Y[2] - array_Y[0]), fabs(array_Y[3] - array_Y[1]));
}

double Calculation_new_F(double h, double toll, double err)
{
    if (fabs(err) > pow(10.0, -18.0))
    {
        return h * min(1.5, max(0.7, 0.98 * pow(toll / err, 1.0 / 6)));
    }
    return h;
}

double eigenvalues(double x, double y)
{
    double lambda_1 = 0;
    double lambda_2 = 1 - (4.0 / 3) * x * (1 + x);
    return fabs(lambda_1) > fabs(lambda_2) ? lambda_1 : lambda_2;
}

void runge_kutta_adaptive(double beginx, double beginz, double h, double toll, std::string file_data, std::string file_delta)
{
    double next_x1 = beginx, next_z1 = beginz, next_x2 = beginx, next_z2 = beginz;
    double err = 0.0, cur_h = 0.0, delta = 0.0;
    double array_Y[4];
    double Value_y[4];
    double t = 0.0, T = 0.0;
    int index = 0;

    std::ofstream out, out_delta;
    out.open(file_data);
    out_delta.open(file_delta);

    out << "t;x;z\n";
    out_delta << "delta\n";

    out << t << ";" << next_x1 << ";" << next_z1 << "\n";
    out_delta << delta << "\n";

    while (t < 10 * M_PI)
    {
        runge_kutta(next_x1, next_z1, h, t, array_Y);
        err = def_err(array_Y);

        do
        {
            h = Calculation_new_F(h, toll, err);
            runge_kutta(next_x1, next_z1, h, t, array_Y);
            err = def_err(array_Y);
        } while (err > toll);

        next_x2 = array_Y[0];
        next_z2 = array_Y[1];

        if (next_z1 * next_z2 < 0)
        {
            if (next_x1 > next_x2)
            {
                t += h;
            }
        }
        else
        {
            t += h;
        }

        delta = err + delta * exp(h * eigenvalues(next_x1, next_z1));
        out_delta << delta << "\n";
        next_x1 = next_x2;
        next_z1 = next_z2;
        out << t << ";" << next_x2 << ";" << next_z2 << "\n";
    }

    std::cout << delta << "\n";

    out.close();
    out_delta.close();
}

int main()
{
    runge_kutta_adaptive(1.0, 0.0, 0.1, 1e-7, "task8-t7.csv", "task8-t7-delta.csv");
    runge_kutta_adaptive(1.0, 0.0, 0.1, 1e-8, "task8-t8.csv", "task8-t8-delta.csv");
    runge_kutta_adaptive(1.0, 0.0, 0.1, 1e-9, "task8-t9.csv", "task8-t9-delta.csv");
    return 0;
}