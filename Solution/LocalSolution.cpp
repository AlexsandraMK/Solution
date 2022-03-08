#include "LocalSolution.h"

double calc_f(Knot coord, double time)   // Вычисление f
{
    //return 0.; // Не зависит от пространства и времени
    //return 0.; // Не зависит от времени
    //return 1.; // Не зависит от пространства
    //return 1.; // Зависит и от пространства и от времени
    //return 2. + 2. * time; // Зависит от времени квадратично
    //return 6 * time + 3. * time * time; // Зависит от времени кубически
    return 12. * time * time + 4. * time * time * time;
}



double prime_by_var(int what, double* Knot, double ksi, double etta, double tetha)
{
    double x = 0.;
    for (int i = 1; i <= 8; i++)
    {
        x += Knot[i - 1] * d_phi(i, what, ksi, etta, tetha);
    }
    return x;
}







double* calc_grad(int index, double ksi, double etta, double tetha)
{
    return new double[3]{ d_phi(index,1,ksi,etta,tetha), d_phi(index,2,ksi,etta,tetha), d_phi(index,3,ksi,etta,tetha) };
}



