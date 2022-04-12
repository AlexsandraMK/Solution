#include "IKE.h"
#include "SolverSLAE.h"



Hexagon::Hexagon()
{
    globalNumsKnots.resize(COUNT_KNOTS);
    knots.resize(COUNT_KNOTS);
    lambda = sigma = hi = 0;
    iterKnots = 0;
}

int Hexagon::GetCountKnots()
{
    return COUNT_KNOTS;
}

vector<vector<double>> Hexagon::CalcLocalM()
{
    vector<vector<double>> M;
    M.resize(COUNT_KNOTS);

    for (int i = 0; i < COUNT_KNOTS; i++)
    {
        M[i].resize(COUNT_KNOTS);
        for (int j = 0; j < COUNT_KNOTS; j++)
            M[i][j] = Integrate(Mij, i, j);
    }
        

    return M;
}

vector<vector<double>> Hexagon::CalcLocalG()
{
    vector<vector<double>> G;
    G.resize(COUNT_KNOTS);

    for (int i = 0; i < COUNT_KNOTS; i++)
    {
        G[i].resize(COUNT_KNOTS);
        for (int j = 0; j < COUNT_KNOTS; j++)
            G[i][j] = Integrate(Gij, i, j);
    }


    return G;
}


int Hexagon::CalcMu(int ind)    {   return (ind % 2) + 1;   }

int Hexagon::CalcNu(int ind)    {   return (ind / 2) % 2 + 1;   }

int Hexagon::CalcZeta(int ind) {   return ind / 4 + 1;   }

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
    double diffPhi = 0;
    int indW = 0;

    switch (axis)
    {
    case ksi:
        indW = CalcMu(ind);
        diffPhi = CalcW(CalcNu(ind), integrationVars[etta]) * CalcW(CalcZeta(ind), integrationVars[theta]);
        break;
    case etta:
        indW = CalcNu(ind);
        diffPhi = CalcW(CalcMu(ind), integrationVars[ksi])  * CalcW(CalcZeta(ind), integrationVars[theta]);
        break;
    case theta:
        indW = CalcZeta(ind);
        diffPhi = CalcW(CalcMu(ind), integrationVars[ksi]) * CalcW(CalcNu(ind), integrationVars[etta]);
        break;
    }

    diffPhi *= indW == 1 ? -1. : 1.;

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
        diffPhiByKsi[i] = DifferentiationPhi(i, ksi, integrationVars);
        diffPhiByEtta[i] = DifferentiationPhi(i, etta, integrationVars);
        diffPhiByTheta[i] = DifferentiationPhi(i, theta, integrationVars);
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

    return grad;
}


double Hexagon::SolveInPoint(Knot knot, vector<double> q)
{
    vector<double> integrationVar = {0.5,0.5,0.5};


    double eps = 1e-16;
    double nkvz = 1e-15;
    
    while (nkvz >= eps)
    {
        
        vector<vector<double>> Jacobian = CalcJacobian(integrationVar);

        double temp;
        for(int i = 0; i < Jacobian.size(); i++)
            for (int j = i; j < Jacobian.size(); j++)
            {
                temp = Jacobian[i][j];
                Jacobian[i][j] = Jacobian[j][i];
                Jacobian[j][i] = temp;
            }

        vector<double> F;
        F.resize(Jacobian.size(),0);

        for (int i = 0; i < COUNT_KNOTS; i++)
        {
            F[0] -= knots[i].x * CalcPhi(i, integrationVar);
            F[1] -= knots[i].y * CalcPhi(i, integrationVar);
            F[2] -= knots[i].z * CalcPhi(i, integrationVar);
        }

        F[0] += knot.x;
        F[1] += knot.y;
        F[2] += knot.z;

        GaussSolverSLAE* solver = new GaussSolverSLAE(Jacobian, F);
        solver->Solve();

        for (int i = 0; i < integrationVar.size(); i++) integrationVar[i] += solver->x[i];
 
        F.resize(Jacobian.size(), 0);
        for (int i = 0; i < COUNT_KNOTS; i++)
        {
            F[0] -= knots[i].x * CalcPhi(i, integrationVar);
            F[1] -= knots[i].y * CalcPhi(i, integrationVar);
            F[2] -= knots[i].z * CalcPhi(i, integrationVar);
        }

        F[0] += knot.x;
        F[1] += knot.y;
        F[2] += knot.z;
        //vector<double> F1 = MultMatrByVect(Jacobian, integrationVar);
        //vector<double> nvkzVec = { 0, 0 , 0 };
        nkvz = 0;
        for (int i = 0; i < F.size(); i++) nkvz += fabs(F[i]);
        //nkvz = sqrt(CalcScalar(nvkzVec, nvkzVec)) / sqrt(CalcScalar(F,F));

    }
    
    vector<double> phiVec = {
        CalcPhi(0, integrationVar),
        CalcPhi(1, integrationVar),
        CalcPhi(2, integrationVar),
        CalcPhi(3, integrationVar),
        CalcPhi(4, integrationVar),
        CalcPhi(5, integrationVar),
        CalcPhi(6, integrationVar),
        CalcPhi(7, integrationVar),
    };


    return q[globalNumsKnots[0]] * phiVec[0] +
        q[globalNumsKnots[1]] * phiVec[1] +
        q[globalNumsKnots[2]] * phiVec[2] +
        q[globalNumsKnots[3]] * phiVec[3] +
        q[globalNumsKnots[4]] * phiVec[4] +
        q[globalNumsKnots[5]] * phiVec[5] +
        q[globalNumsKnots[6]] * phiVec[6] +
        q[globalNumsKnots[7]] * phiVec[7];


}

bool Hexagon::IsIn(Knot knot)
{
    Triangle* leftTrBase = new Triangle();
    leftTrBase->SetGlobalKnotNum(0, knots[0]);
    leftTrBase->SetGlobalKnotNum(2, knots[2]);
    leftTrBase->SetGlobalKnotNum(1, knots[1]);

    Triangle* rightTrBase = new Triangle();
    rightTrBase->SetGlobalKnotNum(1, knots[1]);
    rightTrBase->SetGlobalKnotNum(2, knots[2]);
    rightTrBase->SetGlobalKnotNum(3, knots[3]);

    if((leftTrBase->IsIn(knot) || rightTrBase->IsIn(knot)) &&
        knots[0].z <= knot.z && knot.z <= knots[COUNT_KNOTS - 1].z) return true;

    return false;
}


