#include "SLAU.h"
#include <iostream>
#include <fstream>
#include "InputFuncs.h"
#include <set>


SLAU::SLAU(InitialData* data)
{
    int slauSize = data->knots.size();
    knots.resize(data->knots.size());
    for (int i = 0; i < data->knots.size(); i++)
        knots[i] = data->knots[i];


    A.resize(slauSize);
    for (int i = 0; i < slauSize; i++)
        A[i].resize(slauSize);
    q.resize(slauSize);
    u.resize(slauSize);
    d.resize(slauSize);

    M = AssemblingGlobalM(data);
    /*WriteMatrix(M);*/
    G = AssemblingGlobalG(data);
   /* WriteMatrix(G);*/
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
    for (int iKE = 0; iKE < nKEs; iKE++)
    {
        IKE* ke = data->KEs[iKE];
        int countKnots = ke->GetCountKnots();
        for (int j = 0; j < countKnots; j++)
        {
            int globalJ = ke->globalNumsKnots[j];
            for (int k = 0; k <countKnots; k++)
            {
                int globalK = ke->globalNumsKnots[k];
                // эллиптическая задача
                /*A[globalJ][globalK] = ke->hi * M[globalJ][globalK] + ke->lambda * G[globalJ][globalK];*/
                A[globalJ][globalK] =
                    ke->hi * M[globalJ][globalK] *
                    2 * ((timeToCalc[3] - timeToCalc[2]) + (timeToCalc[3] - timeToCalc[1]) + (timeToCalc[3] - timeToCalc[0]))
                    / ((timeToCalc[3] - timeToCalc[0]) * (timeToCalc[3] - timeToCalc[1]) * (timeToCalc[3] - timeToCalc[2])) +
                    ke->sigma * M[globalJ][globalK]
                    * ((timeToCalc[3] - timeToCalc[1]) * (timeToCalc[3] - timeToCalc[2])
                    +  (timeToCalc[3] - timeToCalc[0]) * (timeToCalc[3] - timeToCalc[2])
                    +  (timeToCalc[3] - timeToCalc[0]) * (timeToCalc[3] - timeToCalc[1]))
                    / ((timeToCalc[3] - timeToCalc[0]) * (timeToCalc[3] - timeToCalc[1]) * (timeToCalc[3] - timeToCalc[2])) +
                    ke->lambda * G[globalJ][globalK];
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

    vector<double>Mq_j3, Mq_j2, Mq_j1;

    Mq_j3 = MultMatrByVect(M, scheme->q[0]);
    Mq_j2 = MultMatrByVect(M, scheme->q[1]);
    Mq_j1 = MultMatrByVect(M, scheme->q[2]);

    int nKEs = data->KEs.size();
    for (int iKE = 0; iKE < nKEs; iKE++)
    {
        IKE* ke = data->KEs[iKE];
        int countKnots = ke->GetCountKnots();
        for (int j = 0; j < countKnots; j++)
        {
            int globalJ = ke->globalNumsKnots[j];
            /*d[globalJ] = b[globalJ];*/
            d[globalJ] = b[globalJ]
                - ke->sigma * Mq_j3[globalJ] * ((timeToCalc[3] - timeToCalc[1]) * (timeToCalc[3] - timeToCalc[2])) /
                                           ((timeToCalc[0] - timeToCalc[1]) * (timeToCalc[0] - timeToCalc[2]) * (timeToCalc[0] - timeToCalc[3]))
                - ke->sigma * Mq_j2[globalJ] * ((timeToCalc[3] - timeToCalc[0]) * (timeToCalc[3] - timeToCalc[2])) /
                                           ((timeToCalc[1] - timeToCalc[0]) * (timeToCalc[1] - timeToCalc[2]) * (timeToCalc[1] - timeToCalc[3]))
                - ke->sigma * Mq_j1[globalJ] * ((timeToCalc[3] - timeToCalc[0]) * (timeToCalc[3] - timeToCalc[1])) /
                                                ((timeToCalc[2] - timeToCalc[0]) * (timeToCalc[2] - timeToCalc[1]) * (timeToCalc[2] - timeToCalc[3]))
                - ke->hi * Mq_j3[globalJ] * 2 * ((timeToCalc[3] - timeToCalc[1]) + (timeToCalc[3] - timeToCalc[2])) /
                                                ((timeToCalc[0] - timeToCalc[3]) * (timeToCalc[0] - timeToCalc[1]) * (timeToCalc[0] - timeToCalc[2]))
                - ke->hi * Mq_j2[globalJ] * 2 * ((timeToCalc[3] - timeToCalc[0]) + (timeToCalc[3] - timeToCalc[2])) /
                                                ((timeToCalc[1] - timeToCalc[0]) * (timeToCalc[1] - timeToCalc[2]) * (timeToCalc[1] - timeToCalc[3]))
                - ke->hi * Mq_j1[globalJ] * 2 * ((timeToCalc[3] - timeToCalc[1]) + (timeToCalc[3] - timeToCalc[0])) /
                                                ((timeToCalc[2] - timeToCalc[0]) * (timeToCalc[2] - timeToCalc[1]) * (timeToCalc[2] - timeToCalc[3]));
        }
    }


}

vector<double> SLAU::AssemblingGlobalF(InitialData* data, double time)
{
    vector<double> f;
    f.resize(data->knots.size());
    int nKEs = data->KEs.size();
    for (int iKE = 0; iKE < nKEs; iKE++)
    {
        IKE* ke = data->KEs[iKE];
        int countKnotsInKe = ke->GetCountKnots();
        vector<double> fInKnots;
        fInKnots.resize(ke->GetCountKnots());
        for (int j = 0; j < countKnotsInKe; j++)
        {
            int globalJ = ke->globalNumsKnots[j];
            fInKnots[j] = GetF(data->knots[globalJ], time);
        }

        vector<double> localF = ke->CalcLocalF(fInKnots);
        for (int j = 0; j < countKnotsInKe; j++)
        {
            int globalJ = ke->globalNumsKnots[j];
            f[globalJ] += localF[j];
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
            int global_num_coord = data->bounds[i].globalNum[j];
            //for (int k = 0; k < data->knots.size(); k++)
            //{
            //    A[global_num_coord][k] = 0;
            //}

            //A[global_num_coord][global_num_coord] = 1.;
            //d[global_num_coord] = u[global_num_coord];
            ////d[global_num_coord] = 0.;     // Для решения


            // БОЛШИМ ЧИСЛОМ
            A[global_num_coord][global_num_coord] = 10e+14;
            d[global_num_coord] = u[global_num_coord] * 10e+14;
            //d[global_num_coord] = 0.* 10e+30;     // Для решения
        }
    }
}

void SLAU::SolveSLAU(InitialData* data, TimeScheme* scheme)
{
    u.resize(data->knots.size(), 0.);
    d.resize(data->knots.size(), 0.);
    for (int i = 0; i < data->knots.size(); i++) A[i].resize(data->knots.size(), 0.);


    CalcU(data, scheme->time[scheme->time.size() - 1]);
    CalcA(data, scheme);
    CalcD(data, scheme);
    CalcFirstBoundaryConditions(data, scheme->time[scheme->time.size() - 1]);
    //WriteMatrix(A);
    LOC();
}

void SLAU::CalcU(InitialData* data, double time)
{
    int slauSize = u.size();
    for (int j = 0; j < slauSize; j++)
    {
        u[j] = GetU(data->knots[j], time);
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
        out << i + 1 << "| ";
        out.width(15);
        out << knots[i].x << "| ";
        out.width(15);
        out << knots[i].y << "| ";
        out.width(15);
        out << knots[i].z << "| ";
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
                    if (diffResult > 1e-14)
                    {
                        out.setf(ios::left);
                        out.width(15);
                        out << knot->x << " ";
                        out.width(15);
                        out << knot->y << " ";
                        out.width(15);
                        out << knot->z << " ";

                        out.width(15);
                        double result = data->KEs[iKe]->SolveInPoint(*knot, q);
                        out << result << " ";

                        out.width(15);
                        double trueResult = GetU(*knot, time);
                        out << trueResult;

                        out.width(15);
                        double diffResult = fabs(trueResult - result);
                        out << diffResult;

                        out << endl;
                    }
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
    z.clear();

    for (int i = 1; i <= 25; i++)
    {
        x.insert(xbeg + i * hx / 25);
        y.insert(ybeg + i * hy / 25);
        z.insert(2.5);
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

    ofstream out(str, ios_base::out | ios_base::app);

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
    out.close();

    //for (int i = 0; i < areaKnots.size(); i++)
    //{
    int i = 0;
    for (set<double> ::iterator iz = z.begin(); iz != z.end(); iz++)
    {
        for (set<double> ::iterator ix = x.begin(); ix != x.end(); ix++)
        {
            for (set<double> ::iterator iy = y.begin(); iy != y.end(); iy++)
            {
                Knot* knot = new Knot(*ix, *iy, *iz);
                //int iKe = FindIKe(data, areaKnots[i]);
                int iKe = FindIKe(data, knot);

                if (iKe >= 0)
                {
                    out.open(str, ios_base::out | ios_base::app);
                    out.setf(ios::left);
                    //out << "| ";
                    //out.width(15);
                    //out << i + 1 << "| ";
                    out.width(15);
                    out << knot->x << " ";
                    out.width(15);
                    out << knot->y << " ";
                    //out.width(15);
                    //out << knot->z << "| ";
                    //out.width(15);

                    double result = data->KEs[iKe]->SolveInPoint(*knot, q);
                    if (fabs(result) < 1e-4) result = 0;
                    out << result;

                    //if (areaKnots[i]->x == 1) out << " " << data->KEs[16]->SolveInPoint(*areaKnots[i], q);
                    //out << "| ";
                    out << endl;
                    delete knot;
                    out.close();
                }
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




