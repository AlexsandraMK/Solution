#include "SolverSLAE.h"

GaussSolverSLAE::GaussSolverSLAE(std::vector<std::vector<double>> A_, std::vector<double> F_)
{
    A = A_;
    F = F_;
    x.resize(F.size());
}

void GaussSolverSLAE::Solve()
{
    double temp1 = 0, temp2 = 0;
    for (int i = 0; i < A.size(); i++)
    {
        int j = i;
        for (; j < A.size() && !A[j][i]; j++);
        if (j == A.size()) break;

        if (j != i)
        {
            SwapLinesInMatrix(A, i, j);
            SwapLinesInVector(F, i, j);
        }

        
        F[i] /= A[i][i];
        for (int k = A.size() - 1; k >= 0; k--) A[i][k] /= A[i][i];

        for (int k = i + 1; k < A.size(); k++)
        {
            double a = A[k][i];
            F[k] -= F[i] * a;
            for (int j = 0; j < A.size(); j++)  A[k][j] -= A[i][j] * a;
        }
    }

    x[2] = F[2] / A[2][2];
    x[1] = (F[1] - x[2] * A[1][2]) / A[1][1];
    x[0] = (F[0] - x[2] * A[0][2] - x[1] * A[0][1]) / A[0][0];

}

