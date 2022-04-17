#include "IntegrationFunc.h"

double Integrate(std::function<double(double, double, double, int, int, double*, double*, double*)> f, int i, int j, double* x, double* y, double* z)
{
    const int nKnot = 5; // ���������� ����� ������

    double xj[nKnot] = { -sqrt(5. + 2. * (sqrt(10. / 7.))) / 3., -sqrt(5. - 2. * (sqrt(10. / 7.))) / 3.,	// �������� ��������� �� �����
               0. , sqrt(5. - 2. * (sqrt(10. / 7.))) / 3. , sqrt(5. + 2. * (sqrt(10. / 7.))) / 3. };

    double qj[nKnot] = { (322. - 13. * sqrt(70.)) / 900., (322. + 13. * sqrt(70.)) / 900., 128. / 225.,	// ���� �����
                   (322. + 13. * sqrt(70.)) / 900., (322. - 13. * sqrt(70.)) / 900. };
    // ����
    double hX = 1. / 2.,
        hY = 1. / 2.,
        hZ = 1. / 2.,
        // ������
        cX = 1. / 2.,
        cY = 1. / 2.,
        cZ = 1. / 2.;

    double result = 0.;
    for (int ix = 0; ix < nKnot; ix++)
        for (int iy = 0; iy < nKnot; iy++)
            for (int iz = 0; iz < nKnot; iz++)
                result += qj[ix] * qj[iy] * qj[iz] * (f(cX + xj[ix] * hX, cY + xj[iy] * hY, cZ + xj[iz] * hZ, i + 1, j + 1, x, y, z));
    return 1. * 1. * 1. * result / 8.; // ���������������
}

double Integrate(std::function<double(std::vector<double> integrationVar, int i, int j)> func, int i, int j)
{
    const int countGaussKnot = 3;//5; // ���������� ����� ������

    // �������� ��������� �� �����
    double xj[countGaussKnot] = { .7745966692414833, 0., -.7745966692414833 };
    /*{ -sqrt(5. + 2. * (sqrt(10. / 7.))) / 3., -sqrt(5. - 2. * (sqrt(10. / 7.))) / 3.,	
    0. , sqrt(5. - 2. * (sqrt(10. / 7.))) / 3. , sqrt(5. + 2. * (sqrt(10. / 7.))) / 3. };*/

    // ���� �����
    double qj[countGaussKnot] = { .55555555555555555, .8888888888888888, .55555555555555555 };
    //{ (322. - 13. * sqrt(70.)) / 900., (322. + 13. * sqrt(70.)) / 900., 128. / 225.,	    
    //(322. + 13. * sqrt(70.)) / 900., (322. - 13. * sqrt(70.)) / 900. };

    std::vector<double> axisSteps;
    std::vector<double> axisCenters;

    // ������ ���� � ������ = 1/2, �.�. �������������� ���������� �� ���� 1x1x1
    axisSteps.resize(3, 1. / 2.);
    axisCenters.resize(3, 1. / 2.);

    double result = 0.;
    std::vector<double> integrationVar;
    for (int ix = 0; ix < countGaussKnot; ix++)
        for (int iy = 0; iy < countGaussKnot; iy++)
            for (int iz = 0; iz < countGaussKnot; iz++)
            {
                integrationVar = { axisCenters[ksi] + xj[ix] * axisSteps[ksi], axisCenters[etta] + xj[iy] * axisSteps[etta], axisCenters[theta] + xj[iz] * axisSteps[theta] };
                result += qj[ix] * qj[iy] * qj[iz] * (func(integrationVar, i, j));
            }

    return result / 8.; // ���������������
}



