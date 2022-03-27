#include <functional>
#include <vector>



double Integrate(std::function<double(double, double, double, int, int, double*, double*, double*)> f, int i, int j, double* x, double* y, double* z);


enum NewAxis {
    ksi = 0,
    etta,
    theta
};

double Integrate(std::function<double(std::vector<double> integrationVar, int i, int j)> func, int i, int j);