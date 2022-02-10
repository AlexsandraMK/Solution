#include <vector>
#include <functional>
using namespace std;

struct Knot // Координаты узлов
{
    double x, y, z;
};

struct LocalArea //локальные области
{
    double h_x, h_y, h_z; // Длины сторон
    double lambda, sigma, hi; // lambda и gamma в области
    int size;
    int* globalNum;  // Номера узлов
    double* x, * y, * z;
};

const int N_TRIANGLE_PRIZM_KNOTS = 6;

const int N_HEXAGON_KNOTS = 8;

void calc_local_M(LocalArea local, double** M);

void calc_local_G(LocalArea local, double** G);

void calc_local_F(vector<Knot> knots, LocalArea local, double F[8], double time);

double Integrate(function<double(double, double, double, int, int, double*, double*, double*)> f, int i, int j, double* x, double* y, double* z);