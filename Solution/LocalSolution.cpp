#include "LocalSolution.h"

double calc_f(Knot coord, double time)   // ���������� f
{
    //return 0.; // �� ������� �� ������������ � �������
    //return 0.; // �� ������� �� �������
    //return 1.; // �� ������� �� ������������
    //return 1.; // ������� � �� ������������ � �� �������
    //return 2. + 2. * time; // ������� �� ������� �����������
    //return 6 * time + 3. * time * time; // ������� �� ������� ���������
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



