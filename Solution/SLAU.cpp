#include "SLAU.h"
#include <iostream>
#include <fstream>
#include "InputFuncs.h"
#include <set>


SLAU::SLAU(InitialData* data)
{
    slauSize = data->knots.size()*3;
    knots.resize(data->knots.size());
    for (int i = 0; i < data->knots.size(); i++)
        knots[i] = data->knots[i];


    A.resize(slauSize);
    for (int i = 0; i < slauSize; i++)
        A[i].resize(slauSize);
    q.resize(slauSize);
    qx.resize(data->knots.size());
    qy.resize(data->knots.size());
    qz.resize(data->knots.size());
    u.resize(slauSize);
    d.resize(slauSize);

    M = AssemblingGlobalM(data);

    

    G_xx = AssemblingGlobalG_aa(data, x, x);
    G_yy = AssemblingGlobalG_aa(data, y, y);
    G_zz = AssemblingGlobalG_aa(data, z, z);
    G_xy = AssemblingGlobalG_aa(data, x, y);
    G_xz = AssemblingGlobalG_aa(data, x, z);
    G_yz = AssemblingGlobalG_aa(data, y, z);

    /*WriteMatrix(M);*/
    /*G = AssemblingGlobalG(data);*/
   /* WriteMatrix(G);*/
}

vector<vector<double>>  SLAU::AssemblingGlobalG_aa(InitialData* data, axis a1, axis a2)
{
    vector<vector<double>> globalMatrix;
    int sizeGlobalMatrix = data->knots.size();
    globalMatrix.resize(sizeGlobalMatrix);
    for (int i = 0; i < sizeGlobalMatrix; i++)
        globalMatrix[i].resize(sizeGlobalMatrix);

    int nKEs = data->KEs.size();
    for (int iKE = 0; iKE < nKEs; iKE++)
    {
        IKE* ke = data->KEs[iKE];
        vector<vector<double>> localMatrix = ke->CalcLocalG();

        int countKnotsInKe = ke->GetCountKnots();

        for (int j = 0; j < countKnotsInKe; j++)
        {
            int globalJ = ke->globalNumsKnots[j];
            for (int k = 0; k < countKnotsInKe; k++)
            {
                int globalK = ke->globalNumsKnots[k];
                globalMatrix[globalJ][globalK] += localMatrix[j][k];
            }
        }
    }

    return globalMatrix;
}


void SLAU::LOC()
{
    int slauSize = q.size();
    int i;
    double nvzk = 0., alfa = 0., beta = 0., skp = 0., eps = 9.999999682655226e-030;

    vector<double> z, r, p, f;
    z.resize(slauSize);
    r.resize(slauSize);
    p.resize(slauSize);
    f.resize(slauSize);

    double lastnvzk;
    for (i = 0; i < slauSize; i++)
        q[i] = 0.;
    f = MultMatrByVect(A, q);

    for (i = 0; i < slauSize; i++)
        z[i] = r[i] = d[i] - f[i];

    p = MultMatrByVect(A, z);
    nvzk = sqrt(CalcScalar(r, r)) / sqrt(CalcScalar(d, d));
    unsigned int k;
    for (k = 1; k < 100000 && nvzk > eps; k++)
    {
        lastnvzk = nvzk;
        skp = CalcScalar(p, p);
        alfa = CalcScalar(p, r) / skp;

        for (i = 0; i < slauSize; i++)
        {
            q[i] += alfa * z[i];
            r[i] -= alfa * p[i];
        }

        f = MultMatrByVect(A, r);
        beta = - CalcScalar(p, f) / skp;

        for (i = 0; i < slauSize; i++)
        {
            z[i] = r[i] + beta * z[i];
            p[i] = f[i] + beta * p[i];
        }

        nvzk = sqrt(CalcScalar(r, r)) / sqrt(CalcScalar(d, d));
    }

}

/// <summary>
/// Вычисление глобальных матриц массы пространственной сетки
/// </summary>
/// <param name="data">Входные данные</param>
vector<vector<double>> SLAU::AssemblingGlobalM(InitialData* data)
{
    vector<vector<double>> globalMatrix;
    int sizeGlobalMatrix = data->knots.size();
    globalMatrix.resize(sizeGlobalMatrix);
    for (int i = 0; i < sizeGlobalMatrix; i++)
        globalMatrix[i].resize(sizeGlobalMatrix);

    int nKEs = data->KEs.size();

    for (int iKE = 0; iKE < nKEs; iKE++)
    {
        IKE* ke = data->KEs[iKE];
        vector<vector<double>> localMatrix = ke->CalcLocalM();
        int countKnotsInKe = ke->GetCountKnots();

        for (int j = 0; j < countKnotsInKe; j++)
        {
            int globalJ = ke->globalNumsKnots[j];
            for (int k = 0; k < countKnotsInKe; k++)
            {
                int globalK = ke->globalNumsKnots[k];
                globalMatrix[globalJ][globalK] += localMatrix[j][k];
            }
        }
    }

    return globalMatrix;
}

/// <summary>
/// Вычисление глобальных матриц жесткости пространственной сетки
/// </summary>
/// <param name="data">Входные данные</param>
vector<vector<double>> SLAU::AssemblingGlobalG(InitialData* data)
{
    vector<vector<double>> globalMatrix;
    int sizeGlobalMatrix = data->knots.size();
    globalMatrix.resize(sizeGlobalMatrix);
    for (int i = 0; i < sizeGlobalMatrix; i++)
        globalMatrix[i].resize(sizeGlobalMatrix);

    int nKEs = data->KEs.size();
    for (int iKE = 0; iKE < nKEs; iKE++)
    {
        IKE* ke = data->KEs[iKE];
        vector<vector<double>> localMatrix = ke->CalcLocalG();

        int countKnotsInKe = ke->GetCountKnots();

        for (int j = 0; j < countKnotsInKe; j++)
        {
            int globalJ = ke->globalNumsKnots[j];
            for (int k = 0; k < countKnotsInKe; k++)
            {
                int globalK = ke->globalNumsKnots[k];
                globalMatrix[globalJ][globalK] += localMatrix[j][k];
            }
        }
    }

    return globalMatrix;
}



void SLAU::CalcA(InitialData* data, TimeScheme* scheme) // Вычисление глобальной A
{
    vector<double> timeToCalc = scheme->time;

    int nKEs = data->KEs.size();
    double density, Yung, Poisson;
    double hi, sigma, lambda;

    for (int iKE = 0; iKE < nKEs; iKE++)
    {
        IKE* ke = data->KEs[iKE];
        density = data->coeffs[ke->iCoeff].density;
        Yung = data->coeffs[ke->iCoeff].Yung;
        Poisson = data->coeffs[ke->iCoeff].Poisson;

        hi = density;
        sigma = Yung / (2. * (1. + Poisson));
        lambda = Poisson * Yung / (1. + Poisson) / (1. - 2. * Poisson);

        int countKnots = ke->GetCountKnots();
        for (int j = 0; j < countKnots; j++)
        {
            int globalJ = ke->globalNumsKnots[j];
            for (int k = 0; k <countKnots; k++)
            {
                int globalK = ke->globalNumsKnots[k] ;

                A[globalJ * 3][globalK * 3] =
                    (lambda + 2* sigma) * G_xx[globalJ][globalK] + 
                    sigma * G_yy[globalJ][globalK] +
                    sigma * G_zz[globalJ][globalK] -
                    hi * M[globalJ][globalK] *
                    2 * ((timeToCalc[3] - timeToCalc[2]) + (timeToCalc[3] - timeToCalc[1]) + (timeToCalc[3] - timeToCalc[0]))
                    / ((timeToCalc[3] - timeToCalc[0]) * (timeToCalc[3] - timeToCalc[1]) * (timeToCalc[3] - timeToCalc[2]));
                
                A[globalJ * 3 +1][globalK * 3 +1] =
                    (lambda + 2 * sigma) * G_yy[globalJ][globalK] +
                    sigma * G_xx[globalJ][globalK] +
                    sigma * G_zz[globalJ][globalK] -
                    hi * M[globalJ][globalK] *
                    2 * ((timeToCalc[3] - timeToCalc[2]) + (timeToCalc[3] - timeToCalc[1]) + (timeToCalc[3] - timeToCalc[0]))
                    / ((timeToCalc[3] - timeToCalc[0]) * (timeToCalc[3] - timeToCalc[1]) * (timeToCalc[3] - timeToCalc[2]));
                
                A[globalJ * 3 +2][globalK * 3 +2] =
                    (lambda + 2 * sigma) * G_zz[globalJ][globalK] +
                    sigma * G_xx[globalJ][globalK] +
                    sigma * G_yy[globalJ][globalK] -
                    hi * M[globalJ][globalK] *
                    2 * ((timeToCalc[3] - timeToCalc[2]) + (timeToCalc[3] - timeToCalc[1]) + (timeToCalc[3] - timeToCalc[0]))
                    / ((timeToCalc[3] - timeToCalc[0]) * (timeToCalc[3] - timeToCalc[1]) * (timeToCalc[3] - timeToCalc[2]));
               
                A[globalJ * 3][globalK * 3 +1] = (lambda + sigma) * G_xy[globalJ][globalK];
                
                A[globalJ * 3][globalK * 3 + 2] = (lambda + sigma) * G_xz[globalJ][globalK];
               
                A[globalJ * 3 + 1][globalK * 3 + 2] = (lambda + sigma) * G_yz[globalJ][globalK];
            }

        }
    }

    //cout << endl;
    //for (int i = 0; i < A.size(); i++)
    //    cout << A[13][i] << "\t";
    //cout << endl;
}

void SLAU::CalcD(InitialData* data, TimeScheme* scheme) // Вычисление глобальной A
{
    vector<double> timeToCalc = scheme->time;


    vector<double> b = AssemblingGlobalF(data, timeToCalc[3]);

    vector<double>Mq_j3x, Mq_j2x, Mq_j1x, Mq_j3y, Mq_j2y, Mq_j1y, Mq_j3z, Mq_j2z, Mq_j1z;

    Mq_j3x = MultMatrByVect(M, scheme->qx[0]);
    Mq_j2x = MultMatrByVect(M, scheme->qx[1]);
    Mq_j1x = MultMatrByVect(M, scheme->qx[2]);
    Mq_j3y = MultMatrByVect(M, scheme->qy[0]);
    Mq_j2y = MultMatrByVect(M, scheme->qy[1]);
    Mq_j1y = MultMatrByVect(M, scheme->qy[2]);
    Mq_j3z = MultMatrByVect(M, scheme->qz[0]);
    Mq_j2z = MultMatrByVect(M, scheme->qz[1]);
    Mq_j1z = MultMatrByVect(M, scheme->qz[2]);
    double density;
    double hi;

    int nKEs = data->KEs.size();
    for (int iKE = 0; iKE < nKEs; iKE++)
    {
        IKE* ke = data->KEs[iKE];

        density = data->coeffs[ke->iCoeff].density;
        hi = density;

        int countKnots = ke->GetCountKnots();
        for (int j = 0; j < countKnots; j++)
        {
            int globalJ = ke->globalNumsKnots[j];
            /*d[globalJ] = b[globalJ];*/
            d[globalJ*3] = b[globalJ*3]
                + hi * Mq_j3x[globalJ] * 2 * ((timeToCalc[3] - timeToCalc[1]) + (timeToCalc[3] - timeToCalc[2])) /
                                                ((timeToCalc[0] - timeToCalc[3]) * (timeToCalc[0] - timeToCalc[1]) * (timeToCalc[0] - timeToCalc[2]))
                + hi * Mq_j2x[globalJ] * 2 * ((timeToCalc[3] - timeToCalc[0]) + (timeToCalc[3] - timeToCalc[2])) /
                                                ((timeToCalc[1] - timeToCalc[0]) * (timeToCalc[1] - timeToCalc[2]) * (timeToCalc[1] - timeToCalc[3]))
                + hi * Mq_j1x[globalJ] * 2 * ((timeToCalc[3] - timeToCalc[1]) + (timeToCalc[3] - timeToCalc[0])) /
                                                ((timeToCalc[2] - timeToCalc[0]) * (timeToCalc[2] - timeToCalc[1]) * (timeToCalc[2] - timeToCalc[3]));

            d[globalJ * 3 + 1] = b[globalJ * 3 + 1]
                + hi * Mq_j3y[globalJ] * 2 * ((timeToCalc[3] - timeToCalc[1]) + (timeToCalc[3] - timeToCalc[2])) /
                ((timeToCalc[0] - timeToCalc[3]) * (timeToCalc[0] - timeToCalc[1]) * (timeToCalc[0] - timeToCalc[2]))
                + hi * Mq_j2y[globalJ] * 2 * ((timeToCalc[3] - timeToCalc[0]) + (timeToCalc[3] - timeToCalc[2])) /
                ((timeToCalc[1] - timeToCalc[0]) * (timeToCalc[1] - timeToCalc[2]) * (timeToCalc[1] - timeToCalc[3]))
                + hi * Mq_j1y[globalJ] * 2 * ((timeToCalc[3] - timeToCalc[1]) + (timeToCalc[3] - timeToCalc[0])) /
                ((timeToCalc[2] - timeToCalc[0]) * (timeToCalc[2] - timeToCalc[1]) * (timeToCalc[2] - timeToCalc[3]));

            d[globalJ * 3 + 2] = b[globalJ * 3 + 2]
                + hi * Mq_j3z[globalJ] * 2 * ((timeToCalc[3] - timeToCalc[1]) + (timeToCalc[3] - timeToCalc[2])) /
                ((timeToCalc[0] - timeToCalc[3]) * (timeToCalc[0] - timeToCalc[1]) * (timeToCalc[0] - timeToCalc[2]))
                + hi * Mq_j2z[globalJ] * 2 * ((timeToCalc[3] - timeToCalc[0]) + (timeToCalc[3] - timeToCalc[2])) /
                ((timeToCalc[1] - timeToCalc[0]) * (timeToCalc[1] - timeToCalc[2]) * (timeToCalc[1] - timeToCalc[3]))
                + hi * Mq_j1z[globalJ] * 2 * ((timeToCalc[3] - timeToCalc[1]) + (timeToCalc[3] - timeToCalc[0])) /
                ((timeToCalc[2] - timeToCalc[0]) * (timeToCalc[2] - timeToCalc[1]) * (timeToCalc[2] - timeToCalc[3]));
        }
    }


}

vector<double> SLAU::AssemblingGlobalF(InitialData* data, double time)
{
    vector<double> f;
    f.resize(data->knots.size()*3);
    int nKEs = data->KEs.size();
    for (int iKE = 0; iKE < nKEs; iKE++)
    {
        IKE* ke = data->KEs[iKE];
        int countKnotsInKe = ke->GetCountKnots();
        vector<double> fInKnotsX, fInKnotsY, fInKnotsZ;
        fInKnotsX.resize(ke->GetCountKnots());
        fInKnotsY.resize(ke->GetCountKnots());
        fInKnotsZ.resize(ke->GetCountKnots());
        for (int j = 0; j < countKnotsInKe; j++)
        {
            int globalJ = ke->globalNumsKnots[j];
            fInKnotsX[j] = GetFx(data->knots[globalJ], time);
            fInKnotsY[j] = GetFy(data->knots[globalJ], time);
            fInKnotsZ[j] = GetFz(data->knots[globalJ], time);
        }

        vector<double> localFx = ke->CalcLocalF(fInKnotsX);
        vector<double> localFy = ke->CalcLocalF(fInKnotsY);
        vector<double> localFz = ke->CalcLocalF(fInKnotsZ);

        for (int j = 0; j < countKnotsInKe; j++)
        {
            int globalJ = ke->globalNumsKnots[j]*3;
            f[globalJ] += localFx[j];
            f[globalJ+1] += localFy[j];
            f[globalJ+2] += localFz[j];
        }
    }

    return f;
}

void SLAU::CalcFirstBoundaryConditions(InitialData* data, double time) 
{
    for (int i = 0; i < data->bounds.size(); i++)
    {
        for (int j = 0; j < 4; j++)
        {
            int global_num_coord = data->bounds[i].globalNum[j]*3;
            for (int k = 0; k < data->knots.size(); k++)
            {
                for(int r1 = 0; r1 <3; r1++)
                    for (int r2 = 0; r2 < 3; r2++)
                        A[global_num_coord + r1][k*3 + r2] = 0;
            }

            A[global_num_coord][global_num_coord] = 1.;
            A[global_num_coord + 1][global_num_coord + 1] = 1.;
            A[global_num_coord + 2][global_num_coord + 2] = 1.;
            d[global_num_coord ] = u[global_num_coord];
            d[global_num_coord +1] = u[global_num_coord+1];
            d[global_num_coord+2] = u[global_num_coord+2];

            //d[global_num_coord] = 0.;     // Для решения


            // БОЛШИМ ЧИСЛОМ
            //A[global_num_coord][global_num_coord] = 10e+14;
            //d[global_num_coord] = u[global_num_coord] * 10e+14;
            //d[global_num_coord] = 0.* 10e+14;     // Для решения
        }
    }
}

void SLAU::SolveSLAU(InitialData* data, TimeScheme* scheme)
{
    u.resize(slauSize, 0.);
    d.resize(slauSize, 0.);
    for (int i = 0; i < slauSize; i++) A[i].resize(slauSize, 0.);


    CalcU(data, scheme->time[scheme->time.size() - 1]);
    CalcA(data, scheme);
    CalcD(data, scheme);
    CalcFirstBoundaryConditions(data, scheme->time[scheme->time.size() - 1]);
    //WriteMatrix(A);
    LOC();

    for (int i = 0; i < data->knots.size(); i++)
    {
        qx[i] = q[i * 3];
        qy[i] = q[i * 3 + 1];
        qz[i] = q[i * 3 + 2];
    }

}

void SLAU::CalcU(InitialData* data, double time)
{
    int slauSize = u.size() / 3;
    for (int j = 0; j < slauSize; j++)
    {
        u[j*3] = GetUx(data->knots[j], time);
        u[j*3+1] = GetUy(data->knots[j], time);
        u[j*3+2] = GetUz(data->knots[j], time);
    }
}

void SLAU::WriteResultForSolution(vector<double> q, double time) //функция вывода в консоль
{
    ofstream out("Result.txt", ios_base::out | ios_base::app);

    out << endl << "ВРЕМЯ: " << time << endl;
    out << endl << "Результат в узлах (веса):" << endl;
    out << " ___________________________________________________________________ " << endl;

    out.setf(ios::left);
    out.width(15);
    out << "| № элемента " << "  | ";
    out.width(15);
    out << "x" << "| ";
    out.width(15);
    out << "y" << "| ";
    out.width(15);
    out << "z" << "| ";
    out.width(15);
    out << "u*" << "| ";
    //out.width(15);
    //out << "u" << "| ";
    //out.width(15);
    //out << "|u-u*|" << "|";
    out << endl;
    out << "|----------------|----------------|";
    out << "----------------|";
    out << "----------------|----------------|";
    out << endl;

    int slauSize = q.size();
    for (int i = 0; i < slauSize; i++)
    {
        out.setf(ios::left);
        out << "| ";
        out.width(15);
        out << i + 1 << "| ";
        out.width(15);
        out << knots[i].x << "| ";
        out.width(15);
        out << knots[i].y << "| ";
        out.width(15);
        out << knots[i].z << "| ";
        out.width(15);
        out << q[i] << "| ";
        //cout.width(15);
        //cout << u[i] << "| ";
        //cout.width(15);
        //cout << fabs(q[i] - u[i]) << "| ";
        out << endl;
    }

    out.close();
}

void SLAU::WriteResultForTest(vector<double> q, double time) //функция вывода в консоль
{
    ofstream out("ResultForTest.txt", ios_base::out | ios_base::app);
    
    out << endl << "ВРЕМЯ: " << time << endl;
    out << "Результат в узлах (веса):" << endl;
    out.setf(ios::left);
    out.width(15);
    out << "| № элемента " << "  | ";
    out.width(15);
    out << "x" << "| ";
    out.width(15);
    out << "y" << "| ";
    out.width(15);
    out << "z" << "| ";
    out.width(15);
    out << "u*" << "| ";
    out.width(15);
    out << "u" << "| ";
    out.width(15);
    out << "|u-u*|" << "|";
    out << endl;

    int slauSize = q.size();
    for (int i = 0; i < slauSize; i++)
    {
        out.setf(ios::left);
        out << "| ";
        out.width(15);
        out << i%3 + 1 << "| ";
        out.width(15);
        out << knots[i%3].x << "| ";
        out.width(15);
        out << knots[i%3].y << "| ";
        out.width(15);
        out << knots[i%3].z << "| ";
        out.width(15);
        out << q[i] << "| ";
        out.width(15);
        out << u[i] << "| ";
        out.width(15);
        out << fabs(q[i] - u[i]) << "| ";
        out << endl;
    }
}

void SLAU::SolveInAreaForTest(InitialData* data, double time) 
{
    set<double> x, y, z;

    for (int i = 0; i < knots.size(); i++)
    {
        x.insert(knots[i].x);
        y.insert(knots[i].y);
        z.insert(knots[i].z);
    }

    double xbeg = *(x.begin());
    double ybeg = *(y.begin());
    double zbeg = *(z.begin());

    double hx = *(--x.end()) - xbeg;
    double hy = *(--y.end()) - ybeg;
    double hz = *(--z.end()) - zbeg;

    x.clear();
    y.clear();
    z.clear();

    int nStepsInArea = 25;

    for (int i = 1; i <= nStepsInArea; i++)
    {
        x.insert(xbeg + i * hx / nStepsInArea);
        y.insert(ybeg + i * hy / nStepsInArea);
        z.insert(zbeg + i * hz / nStepsInArea);
    }

    string str = "ResultArea" + std::to_string(time) + ".txt";

    ofstream out(str);
    int i = 0;
    for (set<double> ::iterator iz = z.begin(); iz != z.end(); iz++)
    {
        for (set<double> ::iterator ix = x.begin(); ix != x.end(); ix++)
        {
            for (set<double> ::iterator iy = y.begin(); iy != y.end(); iy++, i++)
            {
                Knot* knot = new Knot(*ix, *iy, *iz);
                int iKe = FindIKe(data, knot);

                if (iKe >= 0)
                {
                    double result = data->KEs[iKe]->SolveInPoint(*knot, q);
                    double trueResult = GetU(*knot, time);
                    double diffResult = fabs(trueResult - result);
                    //if (diffResult > 1e-8)
                    //{
                        out.setf(ios::left);
                        out.width(15);
                        out << knot->x << " ";
                        out.width(15);
                        out << knot->y << " ";
                        out.width(15);
                        out << knot->z << " ";

                        out.width(15);
                        
                        out << result << " ";

                        out.width(15);
                        
                        out << trueResult;

                        out.width(15);
                        
                        out << diffResult;

                        out << endl;
                    //}
                }

                delete knot;
            }
        }

    }

    if (out.is_open())  out.close();
}


void SLAU::SolveInArea( InitialData* data, double time) //функция вывода в консоль
{
    set<double> x, y, z;

    for (int i = 0; i < knots.size(); i++)
    {
        x.insert(knots[i].x);
        y.insert(knots[i].y);
        z.insert(knots[i].z);
    }

    double hx = *(--x.end()) - *(x.begin());
    double hy = *(--y.end()) - *(y.begin());
    double hz = *(--z.end()) - *(z.begin());
    double xbeg = *(x.begin());
    double ybeg = *(y.begin());
    double zbeg = *(z.begin());

    //x.clear();
    //y.clear();
    y.clear();

    for (int i = 1; i <= 25; i++)
    {
        z.insert(zbeg + i * hz / 25);
        x.insert(xbeg + i * hx / 25);
        y.insert(0.);
    }

    /*vector<Knot*> areaKnots;
    for (set<double> ::iterator iz = z.begin(); iz != z.end(); iz++ )
    {
        for (set<double> ::iterator ix = x.begin(); ix != x.end(); ix++ )
        {
            for (set<double> ::iterator iy = y.begin(); iy != y.end(); iy++)
            {
                Knot* knot = new Knot(*ix, *iy, *iz);
                areaKnots.push_back(knot);
            }
        }
    }*/

    string str = "ResultArea"+std::to_string(time) + ".txt";

    ofstream out(str);

    /*out << endl << "ВРЕМЯ: " << time << endl;
    out << endl << "Результат в узлах (веса):" << endl;
    out << " _____________________________________________________________________________________ " << endl;

    out.setf(ios::left);
    out.width(15);
    out << "| № элемента " << "  | ";
    out.width(15);
    out << "x" << "| ";
    out.width(15);
    out << "y" << "| ";
    out.width(15);
    out << "z" << "| ";
    out.width(15);
    out << "u*" << "| ";
    out.width(15);
    out << endl;
    out << "|----------------|";
    out << "----------------|----------------|----------------|";
    out << "----------------|";
    out << endl;*/
    /*out.close();*/

    //for (int i = 0; i < areaKnots.size(); i++)
    //{
    int i = 0;
    for (set<double> ::iterator iy = y.begin(); iy != y.end(); iy++)
    {
        for (set<double> ::iterator ix = x.begin(); ix != x.end(); ix++)
        {
            for (set<double> ::iterator iz = z.begin(); iz != z.end(); iz++)
            {
                Knot* knot = new Knot(*ix, *iy, *iz);
                //int iKe = FindIKe(data, areaKnots[i]);
                int iKe = FindIKe(data, knot);


                /*out.open(str, ios_base::out | ios_base::app);*/
                out.setf(ios::left);
                //out << "| ";
                //out.width(15);
                //out << i + 1 << "| ";
                out.width(15);
                out << knot->x << " ";
                out.width(15);
                out << knot->z << " ";
                //out.width(15);
                //out << knot->z << "| ";
                //out.width(15);

                double result;
                if (iKe < 0) result = 0.;
                else  result = data->KEs[iKe]->SolveInPoint(*knot, q);
                //if (fabs(result) < 9e-5) result = 0;
                out << result;

                //if (areaKnots[i]->x == 1) out << " " << data->KEs[16]->SolveInPoint(*areaKnots[i], q);
                //out << "| ";
                out << endl;
                delete knot;
                /*out.close()*/;
            }
        }

    }
    
    if(out.is_open())out.close();
}

int SLAU::FindIKe(InitialData* data, Knot* knot)
{
    int iKe = 0;
    for (; iKe < data->KEs.size(); iKe++)
    {
        if (data->KEs[iKe]->IsIn(*knot))
        {
            return iKe;
        }
    }

    cout << "Точка - " << knot->x << " " << knot->y << " " << knot->z << " - находится вне расчетной области\n";

    return -1;
}

double SLAU::SolveInPoint(InitialData* data, Knot knot) 
{
    int iKe = 0;
    for (; iKe < data->KEs.size(); iKe++)
    {
        if (data->KEs[iKe]->IsIn(knot))
        {
            return data->KEs[iKe]->SolveInPoint(knot, q);
        }
    }

    cout << "Точка - " << knot.x << " " << knot.y << " " << knot.z << " - находится вне расчетной области";
    return 0;
}




