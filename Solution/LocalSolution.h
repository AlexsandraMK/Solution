#include <vector>
#include <functional>
using namespace std;

struct Knot // ���������� �����
{
    double x, y, z;
};

struct LocalArea //��������� �������
{
    double h_x, h_y, h_z; // ����� ������
    double lambda, sigma, hi; // lambda � gamma � �������
    int size;
    int* globalNum;  // ������ �����
    double* x, * y, * z;
};

const int N_TRIANGLE_PRIZM_KNOTS = 6;

const int N_HEXAGON_KNOTS = 8;

void calc_local_M(LocalArea local, double** M);

void calc_local_G(LocalArea local, double** G);

void calc_local_F(vector<Knot> knots, LocalArea local, double F[8], double time);

double Integrate(function<double(double, double, double, int, int, double*, double*, double*)> f, int i, int j, double* x, double* y, double* z);