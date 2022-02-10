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

void calc_local_M(LocalArea local, double** M) // Вычисление локальной матрицы массы
{
    //double M_1[8][8] =
    //{
    //    {8, 4, 4, 2, 4, 2, 2, 1},
    //    {4, 8, 2, 4, 2, 4, 1, 2},
    //    {4, 2, 8, 4, 2, 1, 4, 2},
    //    {2, 4, 4, 8, 1, 2, 2, 4},
    //    {4, 2, 2, 1, 8, 4, 4, 2},
    //    {2, 4, 1, 2, 4, 8, 2, 4},
    //    {2, 1, 4, 2, 4, 2, 8, 4},
    //    {1, 2, 2, 4, 2, 4, 4, 8}
    //};

    switch (local.size)
    {
    case (8):
        for (int i = 0; i < 8; i++)
            for (int j = 0; j < 8; j++)
                M[i][j] = Integrate(Mij, i, j, local.x, local.y, local.z);
        break;
    case (6):
        break;
    }


}

void calc_local_G(LocalArea local, double** G) // Вычисление локальной матрицы жесткости
{
    //double G_x[8][8] =  // матрица жесткости(x)
    //{
    //    {4, -4, 2, -2, 2, -2, 1, -1},
    //    {-4, 4, -2, 2, -2, 2, -1, 1},
    //    {2, -2, 4, -4, 1, -1, 2, -2},
    //    {-2, 2, -4, 4, -1, 1, -2, 2},
    //    {2, -2, 1, -1, 4, -4, 2, -2},
    //    {-2, 2, -1, 1, -4, 4, -2, 2},
    //    {1, -1, 2, -2, 2, -2, 4, -4},
    //    {-1, 1, -2, 2, -2, 2, -4, 4}
    //};

    //double G_y[8][8] =  // матрица жесткости(y)
    //{
    //    {4, 2, -4, -2, 2, 1, -2, -1},
    //    {2, 4, -2, -4, 1, 2, -1, -2},
    //    {-4, -2, 4, 2, -2, -1, 2, 1},
    //    {-2, -4, 2, 4, -1, -2, 1, 2},
    //    {2, 1, -2, -1, 4, 2, -4, -2},
    //    {1, 2, -1, -2, 2, 4, -2, -4},
    //    {-2, -1, 2, 1, -4, -2, 4, 2},
    //    {-1, -2, 1, 2, -2, -4, 2, 4}
    //};

    //double G_z[8][8] =  // матрица жесткости(z)
    //{
    //    {4, 2, 2, 1, -4, -2, -2, -1},
    //    {2, 4, 1, 2, -2, -4, -1, -2},
    //    {2, 1, 4, 2, -2, -1, -4, -2},
    //    {1, 2, 2, 4, -1, -2, -2, -4},
    //    {-4, -2, -2, -1, 4, 2, 2, 1},
    //    {-2, -4, -1, -2, 2, 4, 1, 2},
    //    {-2, -1, -4, -2, 2, 1, 4, 2},
    //    {-1, -2, -2, -4, 1, 2, 2, 4}
    //};

    for (int i = 0; i < 8; i++)
        for (int j = 0; j < 8; j++)
            G[i][j] = local.lambda * Integrate(Gij, i, j, local.x, local.y, local.z);
    /*local.h_x * local.h_y * local.h_z / 36 *
    (
        1 / (local.h_x * local.h_x) * G_x[indLocalArea][j] +
        1 / (local.h_y * local.h_y) * G_y[indLocalArea][j] +
        1 / (local.h_z * local.h_z) * G_z[indLocalArea][j]
        );*/
}

void mult_matr_by_vect(int size, double** matrix, double* vector, double* result)
{
    for (int i = 0; i < size; i++)
    {
        double sum = 0;
        for (int j = 0; j < size; j++)
            sum += matrix[i][j] * vector[j];
        result[i] = sum;
    }
}

void calc_local_F(vector<Knot> knots, LocalArea local, double F[8], double time) // Вычисление локального вектора правой части
{
    double** M_1 = new double*[8];
    for (int i = 0; i < 8; i++) M_1[i] = new double[8];
    calc_local_M(local, M_1);

    /*   double M_1[8][8] =
       {
           {8, 4, 4, 2, 4, 2, 2, 1},
           {4, 8, 2, 4, 2, 4, 1, 2},
           {4, 2, 8, 4, 2, 1, 4, 2},
           {2, 4, 4, 8, 1, 2, 2, 4},
           {4, 2, 2, 1, 8, 4, 4, 2},
           {2, 4, 1, 2, 4, 8, 2, 4},
           {2, 1, 4, 2, 4, 2, 8, 4},
           {1, 2, 2, 4, 2, 4, 4, 8}
       };*/

    double f[8]{};

    for (int i = 0; i < 8; i++)
        f[i] = calc_f(knots[local.globalNum[i] - 1], time);

    mult_matr_by_vect(8, M_1, f, F);

    for (int i = 0; i < 8; i++)
        F[i] *= local.h_x * local.h_y * local.h_z / 216;

}

int mu(int index)
{
    return ((index - 1) % 2) + 1;
}

int v(int index)
{
    return ((index - 1) / 2) % 2 + 1;
}

int nu(int index)
{
    return (index - 1) / 4 + 1;
}

double W(int index, double alpha)
{
    switch (index)
    {
    case 1: return 1. - alpha;
    case 2: return alpha;
    }
}

double d_phi(int index, int what, double ksi, double etta, double tetha)
{
    double d_phi = 0.;
    switch (what)
    {
    case 1:    // ksi
    {
        d_phi = W(v(index), etta) * W(nu(index), tetha);
        if (mu(index) == 1) d_phi *= -1;
        break;
    }
    case 2:     // etta
    {
        d_phi = W(mu(index), ksi) * W(nu(index), tetha);
        if (v(index) == 1) d_phi *= -1;
        break;
    }
    case 3:     // tetha
    {
        d_phi = W(mu(index), ksi) * W(v(index), etta);
        if (nu(index) == 1) d_phi *= -1;
        break;
    }
    }

    return d_phi;
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

double phi(int index, double ksi, double etta, double tetha)
{
    return  W(mu(index), ksi) * W(v(index), etta) * W(nu(index), tetha);
}

double det_Jacobian(double* x, double* y, double* z, double ksi, double etta, double tetha)
{
    return prime_by_var(1, x, ksi, etta, tetha) * prime_by_var(2, y, ksi, etta, tetha) * prime_by_var(3, z, ksi, etta, tetha)
        + prime_by_var(3, x, ksi, etta, tetha) * prime_by_var(1, y, ksi, etta, tetha) * prime_by_var(2, z, ksi, etta, tetha)
        + prime_by_var(2, x, ksi, etta, tetha) * prime_by_var(3, y, ksi, etta, tetha) * prime_by_var(1, z, ksi, etta, tetha)
        - prime_by_var(3, x, ksi, etta, tetha) * prime_by_var(2, y, ksi, etta, tetha) * prime_by_var(1, z, ksi, etta, tetha)
        - prime_by_var(1, x, ksi, etta, tetha) * prime_by_var(3, y, ksi, etta, tetha) * prime_by_var(2, z, ksi, etta, tetha)
        - prime_by_var(2, x, ksi, etta, tetha) * prime_by_var(1, y, ksi, etta, tetha) * prime_by_var(3, z, ksi, etta, tetha);
}

double Integrate(function<double(double, double, double, int, int, double*, double*, double*)> f, int i, int j, double* x, double* y, double* z)
{
    const int nKnot = 5; // Количество узлов Гаусса

    double xj[nKnot] = { -sqrt(5. + 2. * (sqrt(10. / 7.))) / 3., -sqrt(5. - 2. * (sqrt(10. / 7.))) / 3.,	// Значения координат на узлах
               0. , sqrt(5. - 2. * (sqrt(10. / 7.))) / 3. , sqrt(5. + 2. * (sqrt(10. / 7.))) / 3. };

    double qj[nKnot] = { (322. - 13. * sqrt(70.)) / 900., (322. + 13. * sqrt(70.)) / 900., 128. / 225.,	// Веса узлов
                   (322. + 13. * sqrt(70.)) / 900., (322. - 13. * sqrt(70.)) / 900. };
    // Шаги
    double hX = 1. / 2.,
        hY = 1. / 2.,
        hZ = 1. / 2.,
        // Центры
        cX = 1. / 2.,
        cY = 1. / 2.,
        cZ = 1. / 2.;

    double result = 0.;
    for (int ix = 0; ix < nKnot; ix++)
        for (int iy = 0; iy < nKnot; iy++)
            for (int iz = 0; iz < nKnot; iz++)
                result += qj[ix] * qj[iy] * qj[iz] * (f(cX + xj[ix] * hX, cY + xj[iy] * hY, cZ + xj[iz] * hZ, i + 1, j + 1, x, y, z));
    return 1. * 1. * 1. * result / 8.; // Масштабирование
}

function<double(double, double, double, int, int, double*, double*, double*)> Mij = [](double ksi, double etta, double tetha, int i, int j, double* x, double* y, double* z)
{
    return phi(i, ksi, etta, tetha) * phi(j, ksi, etta, tetha) * det_Jacobian(x, y, z, ksi, etta, tetha);
};

double* calc_grad(int index, double ksi, double etta, double tetha)
{
    return new double[3]{ d_phi(index,1,ksi,etta,tetha), d_phi(index,2,ksi,etta,tetha), d_phi(index,3,ksi,etta,tetha) };
}

function<double(double, double, double, int, int, double*, double*, double*)> Gij = [](double ksi, double etta, double tetha, int i, int j, double* x, double* y, double* z)
{


    double* J_grad_i = new double[3]{};
    double* J_grad_j = new double[3]{};

    double** reversed_Jacobian = new double* [3]
    {
        new double[3]{},
        new double[3]{},
        new double[3]{}
    };


    double** Jacobian = new double* [3]
    {
        new double[3]{prime_by_var(1,x, ksi, etta, tetha),prime_by_var(1,y, ksi, etta, tetha),prime_by_var(1,z, ksi, etta, tetha)},
        new double[3]{prime_by_var(2,x, ksi,etta,tetha) ,prime_by_var(2,y, ksi,etta,tetha),prime_by_var(2,z, ksi,etta,tetha)},
        new double[3]{prime_by_var(3,x, ksi,etta,tetha),prime_by_var(3,y, ksi,etta,tetha),prime_by_var(3,z, ksi,etta,tetha)}
    };

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            double min[4]{};
            int k = 0;
            for (int im = 0; im < 3; im++)
                for (int jm = 0; jm < 3; jm++)
                {
                    if (im != i && jm != j)
                        min[k++] = Jacobian[im][jm];
                }
            reversed_Jacobian[j][i] = pow(-1, i + j + 2) * (min[0] * min[3] - min[1] * min[2]);
        }
    }



    mult_matr_by_vect(3, reversed_Jacobian, calc_grad(i, ksi, etta, tetha), J_grad_i);
    mult_matr_by_vect(3, reversed_Jacobian, calc_grad(j, ksi, etta, tetha), J_grad_j);

    return scalar(3, J_grad_i, J_grad_j) / det_Jacobian(x, y, z, ksi, etta, tetha);
};

