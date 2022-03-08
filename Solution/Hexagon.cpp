#include "IKE.h"



Hexagon::Hexagon()
{
    globalNumsKnots.resize(COUNT_KNOTS);
    lambda = sigma = hi = 0;
    iterKnots = 0;
}

vector<vector<double>> Hexagon::CalcLocalM()
{
    vector<vector<double>> M;
    M.resize(COUNT_KNOTS);

    for (int i = 0; i < 8; i++)
        for (int j = 0; j < 8; j++)
            M[i][j] = Integrate(Mij, i, j);

    return M;
}

vector<vector<double>> Hexagon::CalcLocalG()
{
    vector<vector<double>> G;
    G.resize(COUNT_KNOTS);

    for (int i = 0; i < 8; i++)
        for (int j = 0; j < 8; j++)
            G[i][j] = Integrate(Gij, i, j);

    return G;
}


int Hexagon::CalcMu(int ind)    {   return ((ind - 1) % 2) + 1;   }

int Hexagon::CalcNu(int ind)    {   return ((ind - 1) / 2) % 2 + 1;   }

int Hexagon::CalcZeta(int ind) {   return (ind - 1) / 4 + 1;   }

double Hexagon::CalcW(int ind, double alpha)
{
	switch (ind)
	{
	case 1: return 1. - alpha;
	case 2: return alpha;
	}
}

double Hexagon::CalcPhi(int ind, vector<double> integrationVar)
{
    return  CalcW(CalcMu(ind), integrationVar[ksi]) * CalcW(CalcNu(ind), integrationVar[etta]) * CalcW(CalcZeta(ind), integrationVar[theta]);
}

double Hexagon::DifferentiationPhi(int ind, NewAxis axis, vector<double> integrationVars)
{
    double diffPhi = CalcPhi(ind, integrationVars);
    int indW = 0;

    switch (axis)
    {
    case ksi:
        indW = CalcMu(ind);
        diffPhi /= CalcW(indW, integrationVars[ksi]);
        break;
    case etta:
        indW = CalcNu(ind);
        diffPhi /= CalcW(indW, integrationVars[etta]);
        break;
    case theta:
        indW = CalcZeta(ind);
        diffPhi /= CalcW(indW, integrationVars[theta]);
        break;
    }

    diffPhi *= indW == 1 ? -1 : 1;

    return diffPhi;
}

vector<vector<double>> Hexagon::CalcJacobian(vector<double> integrationVars)
{
    vector<double> x;
    vector<double> y;
    vector<double> z;
    vector<double> diffPhiByKsi;
    vector<double> diffPhiByEtta;
    vector<double> diffPhiByTheta;

    x.resize(COUNT_KNOTS);
    y.resize(COUNT_KNOTS);
    z.resize(COUNT_KNOTS);
    diffPhiByKsi.resize(COUNT_KNOTS);
    diffPhiByEtta.resize(COUNT_KNOTS);
    diffPhiByTheta.resize(COUNT_KNOTS);


    for (int i = 0; i < COUNT_KNOTS; i++)
    {
        x[i] = knots[i].x;
        y[i] = knots[i].y;
        z[i] = knots[i].z;
        diffPhiByKsi[i] = DifferentiationPhi(i + 1, ksi, integrationVars);
        diffPhiByEtta[i] = DifferentiationPhi(i + 1, etta, integrationVars);
        diffPhiByTheta[i] = DifferentiationPhi(i + 1, theta, integrationVars);
    }

    vector<vector<double>> Jacobian =
    {
        {   CalcScalar(x, diffPhiByKsi),    CalcScalar(y, diffPhiByKsi),    CalcScalar(z, diffPhiByKsi)     },
        {   CalcScalar(x, diffPhiByEtta),   CalcScalar(y, diffPhiByEtta),   CalcScalar(z, diffPhiByEtta)    },
        {   CalcScalar(x, diffPhiByTheta),  CalcScalar(y, diffPhiByTheta),  CalcScalar(z, diffPhiByTheta)   }
    };

    return Jacobian;
}

vector<double> Hexagon::CalcGrad(int ind, vector<double> integrationVars)
{
    vector<double> grad = { DifferentiationPhi(ind,ksi, integrationVars),
                            DifferentiationPhi(ind,etta, integrationVars),
                            DifferentiationPhi(ind,theta, integrationVars) };
}


